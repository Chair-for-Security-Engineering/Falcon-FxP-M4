#!/usr/bin/env python3
"""Constant-time scan of a Cortex-M4 image (pqm4 ELF, built with -g3).

Walks the call graph from a root symbol (default crypto_sign_signature) and
checks the instructions whose timing can depend on data:

  hard checks, must be zero on the whole graph:
    UDIV/SDIV (early termination, 2-12 cycles), libgcc arithmetic helpers,
    VFP arithmetic, indirect calls;
  triage, one row per site:
    conditional branches and IT blocks, mapped to their source line and
    classified by the C construct there (loop/if/ternary/logic/switch, or
    'compiler' when the line has none, i.e. the compiler introduced it).
    An IT block is 'alu' when every predicated instruction is a
    single-cycle ALU op, which costs the same whether it executes or not.

A branch on a line with a C condition implements that condition, which the
host TIMECOP run checks (scripts/ct-matrix.sh, also at -O0 where every C
condition becomes a jump), provided test_ct executes the line
(scripts/ct_cover.py checks that). The rows that need a reason here are the
'compiler', 'asm' and 'nosrc' (source not found) ones and the IT blocks
predicating more than ALU ops; an allowlist (lines 'function|file:line|reason')
marks the reviewed ones, any other such row is printed as UNREVIEWED.

  scripts/ct_scan.py ELF [--root SYM] [--allow FILE] [--src DIR] [--all]
"""
import argparse
import collections
import os
import re
import subprocess
import sys

COND = r'(eq|ne|cs|cc|hs|lo|mi|pl|vs|vc|hi|ls|ge|lt|gt|le)'
LIBGCC = re.compile(r'^__(aeabi_|ash|lsh|udiv|div|mod|umod|clz|ctz|mul'
                    r'|add|sub|fix|float|extend|trunc|neg|cmp)')
ALU = re.compile(r'^(mov|mvn|add|adc|sub|sbc|rsb|neg|and|orr|orn|eor|bic|lsl'
                 r'|lsr|asr|ror|cmp|cmn|tst|teq|mul|mla|sxt|uxt|sel)')
# Library code on public data only (sizes, alignment, the board TRNG).
PUBLIC_FUNCS = {'memcpy', 'memset', 'memmove', 'PQCLEAN_randombytes',
                'rng_get_random', 'rng_get_random_blocking', 'randombytes'}


def tool(name):
    return os.environ.get('CROSS', 'arm-none-eabi-') + name


def disassemble(elf):
    out = subprocess.run([tool('objdump'), '-d', '--no-show-raw-insn', elf],
                         capture_output=True, text=True, check=True).stdout
    body = collections.OrderedDict()
    fn = None
    for line in out.splitlines():
        m = re.match(r'^([0-9a-f]+) <(.+)>:$', line)
        if m:
            fn = m.group(2)
            body[fn] = []
            continue
        m = re.match(r'^\s+([0-9a-f]+):\t([^\t]+)\t?(.*)$', line)
        if fn and m:
            body[fn].append((int(m.group(1), 16), m.group(2).strip(),
                             m.group(3).strip()))
    return body


def reachable(body, roots):
    calls = collections.defaultdict(set)
    for f, ins in body.items():
        for _, op, arg in ins:
            t = re.search(r'<([^+>]+)>', arg)
            if op in ('bl', 'b', 'b.n', 'b.w') and t and t.group(1) != f:
                calls[f].add(t.group(1))
    seen, todo = set(), list(roots)
    while todo:
        f = todo.pop()
        if f in seen or f not in body:
            continue
        seen.add(f)
        todo.extend(calls[f])
    return seen


def addr2line(elf, addrs):
    if not addrs:
        return {}
    p = subprocess.run([tool('addr2line'), '-e', elf, '-a'] +
                       [hex(a) for a in addrs],
                       capture_output=True, text=True, check=True)
    lines = p.stdout.splitlines()
    res = {}
    for i in range(0, len(lines) - 1, 2):
        res[int(lines[i], 16)] = lines[i + 1].strip()
    return res


_src_cache = {}


def source_text(loc, srcdirs):
    m = re.match(r'^(.*):(\d+)', loc)
    if not m:
        return None
    path, num = m.group(1), int(m.group(2))
    cands = [path] + [os.path.join(d, os.path.basename(path))
                      for d in srcdirs]
    for c in cands:
        if c not in _src_cache:
            try:
                with open(c, errors='replace') as fh:
                    _src_cache[c] = fh.read().split('\n')
            except OSError:
                _src_cache[c] = None
        if _src_cache[c] is not None and 0 < num <= len(_src_cache[c]):
            return _src_cache[c][num - 1].strip()
    return None


def classify(loc, text):
    if loc.endswith(('.s', '.S')) or re.search(r'\.[sS]:\d+', loc):
        return 'asm'
    if text is None:
        return 'nosrc'
    t = re.sub(r'/\*.*?\*/|//.*$', '', text)
    if re.search(r'\b(for|while)\s*\(|^\s*do\b|\}\s*while\b', t):
        return 'loop'
    if re.search(r'\bswitch\s*\(|\bcase\b', t):
        return 'switch'
    if re.search(r'\bif\s*\(', t):
        return 'if'
    if '?' in t:
        return 'ternary'
    if '&&' in t or '||' in t:
        return 'logic'
    return 'compiler'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('elf')
    ap.add_argument('--root', action='append')
    ap.add_argument('--allow')
    ap.add_argument('--src', action='append', default=[])
    ap.add_argument('--all', action='store_true',
                    help='print every row, not only unreviewed ones')
    a = ap.parse_args()
    roots = a.root or ['crypto_sign_signature']
    body = disassemble(a.elf)
    funcs = reachable(body, roots)

    hard = collections.Counter()
    sites = []   # (func, addr, kind, detail)
    for f in sorted(funcs):
        ins = body[f]
        i = 0
        while i < len(ins):
            addr, op, arg = ins[i]
            t = re.search(r'<([^+>]+)>', arg)
            if op in ('udiv', 'sdiv'):
                hard['div'] += 1
                print(f'HARD div {f} {addr:x} {op} {arg}')
            if op == 'bl' and t and LIBGCC.match(t.group(1)):
                hard['libgcc'] += 1
                print(f'HARD libgcc {f} {addr:x} -> {t.group(1)}')
            if re.match(r'v(add|sub|mul|mla|mls|fma|fms|div|sqrt|cvt|nmul)',
                        op):
                hard['vfp'] += 1
                print(f'HARD vfp {f} {addr:x} {op}')
            if op == 'blx' or (op == 'bx' and arg != 'lr'):
                hard['indirect'] += 1
                print(f'HARD indirect {f} {addr:x} {op} {arg}')
            if f in PUBLIC_FUNCS:
                i += 1
                continue
            if (re.fullmatch(r'b' + COND + r'(\.n|\.w)?', op)
                    or op in ('cbz', 'cbnz')):
                tgt = int(arg.split()[0], 16) if re.match(
                    r'^[0-9a-f]+\b', arg) else None
                back = tgt is not None and tgt <= addr
                sites.append((f, addr, 'br-back' if back else 'br-fwd',
                              f'{op} {arg}'))
            m = re.fullmatch(r'it([te]{0,3})', op)
            if m:
                n = 1 + len(m.group(1))
                pred = [p[1] for p in ins[i + 1:i + 1 + n]]
                alu = all(ALU.match(p) for p in pred)
                sites.append((f, addr, 'it-alu' if alu else 'it-other',
                              f'{op} {arg}: ' + ' '.join(pred)))
            i += 1

    locs = addr2line(a.elf, [s[1] for s in sites])
    allow = {}
    if a.allow:
        with open(a.allow) as fh:
            for line in fh:
                line = line.strip()
                if line and not line.startswith('#'):
                    parts = line.split('|', 2)
                    if len(parts) == 3:
                        allow[(parts[0], parts[1])] = parts[2]
    srcdirs = a.src + [os.path.dirname(os.path.abspath(a.elf))]
    cats = collections.Counter()
    unreviewed = 0
    rows = []
    for f, addr, kind, detail in sites:
        loc = locs.get(addr, '??:0')
        loc = re.sub(r' \(discriminator \d+\)', '', loc)
        short = os.path.basename(loc)
        text = source_text(loc, srcdirs)
        cat = classify(loc, text)
        cats[(kind, cat)] += 1
        reason = allow.get((f, short))
        if reason is None and kind != 'it-alu' and (
                kind == 'it-other' or cat in ('compiler', 'asm', 'nosrc')):
            reason = None
            unreviewed += 1
        elif reason is None:
            reason = 'C condition' if kind != 'it-alu' else 'alu-only'
        rows.append((f, addr, kind, cat, short, text or '', detail, reason))

    print(f'== {a.elf}')
    print(f'roots {",".join(roots)}: {len(funcs)} functions; hard checks: '
          + ' '.join(f'{k}={hard[k]}' for k in
                     ('div', 'libgcc', 'vfp', 'indirect')))
    print('sites by kind/category: ' + ', '.join(
        f'{k}/{c}={n}' for (k, c), n in sorted(cats.items())))
    for f, addr, kind, cat, short, text, detail, reason in rows:
        if a.all or reason is None:
            tag = reason if reason is not None else 'UNREVIEWED'
            print(f'{kind:8s} {cat:8s} {f}|{short} [{addr:x}] {detail}'
                  f'  // {text[:70]}  => {tag}')
    print(f'unreviewed sites: {unreviewed}')
    sys.exit(1 if any(hard.values()) else 0)


if __name__ == '__main__':
    main()
