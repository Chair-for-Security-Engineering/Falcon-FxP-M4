#!/usr/bin/env python3
"""Links the Cortex-M4 scan to the host TIMECOP run: is every C condition
behind a branch of a signing image executed by test_ct?

scripts/ct_scan.py trusts a branch on a line holding a C condition to the
host TIMECOP run (scripts/ct-matrix.sh). That holds only if test_ct runs
the line. For each configuration this builds both pqm4 images, lists their
'C condition' sites (ct_scan.py --all), builds test_ct with --coverage at
-O0, runs it, and reports the sites gcov never saw executed. Each one must
be reviewed and go into scripts/ct_scan.allow. Enters shell-cm4.nix when
the Arm binutils are missing; run from the repository root. Same
configurations as scripts/ct-matrix.sh.

  [CT_COVER_KEYS=20] [CT_COVER_SIGS=20] scripts/ct_cover.py [filter]
"""
import collections
import glob
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

OUT = 'results/ct-cover'
KEYS = os.environ.get('CT_COVER_KEYS', '20')
SIGS = os.environ.get('CT_COVER_SIGS', '20')
# name, pqm4 implementation, pqm4 knobs, host knobs: the paper's board
# configurations (configs.mk) except ieee-base-origshake, plus the fxp split
# row without the 64-bit l10 and ties-up.
PRESETS = [
    ('ieee-base', 'ieee-m4f'), ('ieee-base-fftnew', 'ieee-m4f'),
    ('ieee-split', 'ieee-m4f'), ('ieee-rec', 'ieee-m4f'),
    ('ieee-rec-fftnew', 'ieee-m4f'), ('fxp-split', 'fxp3264-cm4'),
    ('fxp-unfused', 'fxp3264-cm4'), ('fxp-rec', 'fxp3264-cm4'),
    ('fxp-rec-fftnew', 'fxp3264-cm4'), ('fxp-nonsplit', 'fxp3264-cm4'),
    ('fxp-nonsplit-tiesup', 'fxp3264-cm4'),
    ('fxp-nonsplit-tiesup-fftnew', 'fxp3264-cm4'),
]
CONFIGS = [(n, i, [f'CONFIG={n}'], [f'CONFIG={n}']) for n, i in PRESETS] + [
    ('fxp-split-nol10', 'fxp3264-cm4', ['L10_64=0', 'TIES_UP=0'],
     ['BACKEND=fxp32.64', 'L10_64=0', 'TIES_UP=0']),
]


def run(cmd, cwd=None, log=None):
    with open(log, 'w') if log else open(os.devnull, 'w') as fh:
        return subprocess.run(cmd, cwd=cwd, stdout=fh,
                              stderr=subprocess.STDOUT).returncode


def m4_sites(name, impl, knobs):
    """C-condition sites (file, line) of the 512 and 1024 images."""
    sites = collections.defaultdict(set)
    for s in ('512', '1024'):
        tag = f'{name}-{s}'
        if run(['make', '-s', 'pqm4', f'BENCH_IMPL={impl}', f'SCHEME={s}']
               + knobs, log=f'{OUT}/{tag}.m4build.log'):
            print(f'BUILD-FAIL m4 {tag}', flush=True)
            return None
        elf = f'pqm4/elf/crypto_sign_fndsa{s}_{impl}_speed.elf'
        rep = f'{OUT}/{tag}.scan-all.txt'
        run(['python3', 'scripts/ct_scan.py', elf, '--root',
             'crypto_sign_signature', '--src', '.', '--src', 'fpr/fxp32.64',
             '--src', 'fpr/ieee', '--src', 'pqm4/common', '--allow',
             'scripts/ct_scan.allow', '--all'], log=rep)
        for line in open(rep):
            m = re.match(r'^(br-\S+)\s+(\S+)\s+(\S+)\|(\S+):(\d+) .*=> (.*)$',
                         line.rstrip())
            if m and m.group(6) == 'C condition':
                sites[(m.group(4), int(m.group(5)))].add(
                    f'{s}:{m.group(3)}:{m.group(2)}')
    return sites


def host_cover(name, knobs):
    """Lines (file, line) executed by a --coverage host test_ct run."""
    root = os.getcwd()
    work = tempfile.mkdtemp(prefix='ctcov-')
    try:
        subprocess.run(['rsync', '-a', '--exclude', '.git', '--exclude',
                        'pqm4', '--exclude', 'results', '--exclude', '*.o',
                        '--exclude', '*.gcda', '--exclude', '*.gcno',
                        root + '/', work + '/'], check=True)
        extra = '-O0 --coverage -fno-tree-vectorize'
        if 'CONFIG=ieee' in ' '.join(knobs):
            extra += ' -DFNDSA_AVX2=0 -DFNDSA_SSE2=0'
        if run(['make', '-s', 'test_ct', 'CC=gcc', f'EXTRA_CFLAGS={extra}',
                'LDFLAGS=--coverage'] + knobs, cwd=work,
               log=f'{root}/{OUT}/{name}.hostbuild.log'):
            print(f'BUILD-FAIL host {name}', flush=True)
            return None
        if run(['./test_ct', KEYS, SIGS], cwd=work,
               log=f'{root}/{OUT}/{name}.hostrun.log'):
            print(f'RUN-FAIL host {name}', flush=True)
            return None
        done = set()
        for gcda in glob.glob(f'{work}/**/*.gcda', recursive=True):
            p = subprocess.run(['gcov', '-t', '-o', os.path.dirname(gcda),
                                gcda[:-5] + '.c'], cwd=work,
                               capture_output=True, text=True)
            cur = None
            for line in p.stdout.split('\n'):
                m = re.match(r'^\s*-:\s*0:Source:(.*)$', line)
                if m:
                    cur = os.path.basename(m.group(1))
                    continue
                m = re.match(r'^\s*(\d+)\*?:\s*(\d+):', line)
                if m and cur and int(m.group(1)) > 0:
                    done.add((cur, int(m.group(2))))
        return done
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main():
    # The Cortex-M4 shell has the Arm binutils, gcc, gcov and rsync.
    if not shutil.which('arm-none-eabi-objdump') and shutil.which('nix-shell'):
        cmd = ' '.join(shlex.quote(a) for a in [sys.argv[0]] + sys.argv[1:])
        os.execvp('nix-shell', ['nix-shell', '-I',
                               'nixpkgs=' + os.path.abspath('nixpkgs.nix'),
                               'shell-cm4.nix', '--run', cmd])
    os.makedirs(OUT, exist_ok=True)
    only = sys.argv[1] if len(sys.argv) > 1 else ''
    total = missed = 0
    with open(f'{OUT}/summary.txt', 'w') as summ:
        for name, impl, pk, hk in CONFIGS:
            if only and not re.search(only, name):
                continue
            sites = m4_sites(name, impl, pk)
            cov = host_cover(name, hk)
            if sites is None or cov is None:
                summ.write(f'ERROR {name}\n')
                missed += 1
                continue
            miss = sorted(k for k in sites if k not in cov)
            total += len(sites)
            missed += len(miss)
            line = (f'{name}: {len(sites)} C-condition lines in the M4 '
                    f'images, {len(miss)} not executed by test_ct')
            print(line, flush=True)
            summ.write(line + '\n')
            for f, n in miss:
                summ.write(f'    {f}:{n}  {sorted(sites[(f, n)])}\n')
        summ.write(f'TOTAL {total} lines, {missed} not executed\n')
    print(f'TOTAL {total} lines, {missed} not executed', flush=True)
    sys.exit(1 if missed else 0)


if __name__ == '__main__':
    main()
