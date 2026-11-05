/*
 * Precision-evaluation harness, built only for BACKEND=prec. Runs seeded
 * keygen and signing and records the mu and 1/sigma errors from the prec
 * backend's sampler checkpoint (see the usage message for options).
 *
 * Seeds depend only on the seed, logn and run index, so any run can be
 * reproduced alone and the data is the same for any number of shards.
 * With --out, writes PREFIX.stats (exponent histograms and squared-error
 * sums, mergeable with --merge), PREFIX.perkey (per-key sums for the
 * per-key RMSE) and PREFIX.worst (worst case per run); --csv also writes
 * every call to PREFIX.csv.
 *
 * An error 0.xxx * 2^-e is reported as e, so the smallest e is the worst
 * case. The summary gives e statistics, the RMSE, the e histogram, the
 * rounding-divergence counters and, after --merge, the per-key RMSE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <signal.h>
#include <unistd.h>
#include <math.h>

#include "fndsa.h"
#include "knobs.h"
#include "fpr/prec/prec_eval.h"
#include "fpr/prec/target_glue.h"

/* ==================================================================== */

/* Compiled-in variant flags as the result-name token string.
   TWIDDLE_U64 has no token: prec builds always run with it off. */
static const char *
variant_tokens(void)
{
	return
#if FNDSA_SIGN_SPLIT
		"-split"
#elif FNDSA_FXP3264
		"-scalefirst"
#endif
#if FNDSA_SIGN_INVQ_WIDE
		"-invqwide"
#endif
#if FNDSA_SIGN_RECIP_FUSED
		"-recipfused"
#endif
#if FNDSA_SIGN_L10_64
		"-l10div"
#endif
#if FNDSA_SIGN_SQRT_FUSED
		"-sqrtfused"
#endif
#if FNDSA_FXP_TIES_UP
		"-tiesup"
#endif
#if !FNDSA_FXP_CMUL_FUSED
		"-cmulcomposed"
#endif
#if FNDSA_SYMP_ROOT
		"-symproot"
#endif
#if FNDSA_SYMP_TREE
		"-symptree"
#endif
#if FNDSA_FFT_NEW
		"-fftnew"
#endif
#if FNDSA_SIGN_INT_GRAM
		"-intgram"
#endif
#if FNDSA_KGEN_GAMMA
		"-gamma"
#endif
		"";
}

/* The compiled-in flags plus the run-time master choice. */
static const char *
name_tokens(void)
{
	static char buf[256];
	const char *m = getenv("PREC_MASTER");

	snprintf(buf, sizeof buf, "%s%s", variant_tokens(),
		(m != 0 && strcmp(m, "ref") == 0) ? "-refmaster" : "");
	return buf;
}

/* The token string for summaries and sidecars ("none" if no flag). */
static const char *
variant_label(const char *tokens)
{
	return (tokens[0] == 0) ? "none" : tokens;
}

static FILE *csv = 0;      /* per-call rows, opt-in via --csv */
static uint64_t cur_run = 0;
static const char *tree_id = "unknown";

/* Per-run worst case; PREC_EVAL_E_MAX + 1 if every error was exact or
   tiny. */
static int run_min_mu, run_min_sig;

static void
run_worst_reset(void)
{
	run_min_mu = PREC_EVAL_E_MAX + 1;
	run_min_sig = PREC_EVAL_E_MAX + 1;
}

/* Per-call CSV cell: 0 for exact, -1 (row dropped) for huge or tiny,
   as in TWFalcon's precision test. */
static int
csv_code(int e)
{
	if (e == PREC_EVAL_E_NONE) {
		return 0;
	}
	if (e <= 0 || e > PREC_EVAL_E_MAX) {
		return -1;
	}
	return e;
}

static void
eval_sink(void *ctx, int e_mu, int e_sigma)
{
	(void)ctx;
	if (e_mu != PREC_EVAL_E_NONE && e_mu < run_min_mu) {
		run_min_mu = e_mu;
	}
	if (e_sigma != PREC_EVAL_E_NONE && e_sigma < run_min_sig) {
		run_min_sig = e_sigma;
	}
	if (csv != 0) {
		int cm = csv_code(e_mu), cs = csv_code(e_sigma);
		if (cm >= 0 && cs >= 0) {
			fprintf(csv, "%llu,%d,%d\n",
				(unsigned long long)cur_run, cm, cs);
		}
	}
}

static void
alarm_handler(int sig)
{
	(void)sig;
	static const char msg[] =
		"eval_precision: run timed out (broken target arithmetic"
		" spinning in the signing retry loop?)\n";
	ssize_t w = write(2, msg, sizeof msg - 1);
	(void)w;
	_exit(3);
}

/* ==================================================================== */
/* Per-key records. The error bound is an RMSE per key, so each shard
   writes one record per key it touched and the merge folds them.      */

typedef struct {
	uint64_t key, sigs, calls;
	double sq_mu, sq_sig;
	int min_mu, min_sig;
} perkey_rec;

static void
perkey_write(FILE *f, const perkey_rec *k)
{
	fprintf(f, "%llu %llu %llu %a %a %d %d\n",
		(unsigned long long)k->key, (unsigned long long)k->sigs,
		(unsigned long long)k->calls, k->sq_mu, k->sq_sig,
		k->min_mu, k->min_sig);
}

/* ==================================================================== */
/* Summary printing.                                                    */

static const char *
sigma_label(void)
{
	/* sigma is a relative error for every target. */
	return "sigma_rel";
}

static void
print_hist_summary(const char *name, const prec_eval_hist *h,
	uint64_t calls)
{
	uint64_t n = 0, sum_e = 0;
	int min_e = 0, max_e = 0;

	for (int e = 1; e <= PREC_EVAL_E_MAX; e ++) {
		if (h->hist[e] != 0) {
			if (n == 0) {
				min_e = e;
			}
			max_e = e;
			n += h->hist[e];
			sum_e += h->hist[e] * (uint64_t)e;
		}
	}
	printf("%-10s n=%llu exact=%llu huge=%llu tiny=%llu",
		name, (unsigned long long)n, (unsigned long long)h->exact,
		(unsigned long long)h->huge, (unsigned long long)h->tiny);
	if (n != 0) {
		uint64_t half = (n + 1) / 2, acc = 0;
		int med = min_e;
		for (int e = min_e; e <= max_e; e ++) {
			acc += h->hist[e];
			if (acc >= half) {
				med = e;
				break;
			}
		}
		printf("  min(worst)=%d max=%d mean=%.2f median=%d",
			min_e, max_e, (double)sum_e / (double)n, med);
		if (min_e <= 40) {
			printf("  [WARN: at/below the e=40 Falcon"
				" requirement marker]");
		}
	}
	/* Huge errors are not in the histogram, so report them here. */
	if (h->huge != 0) {
		printf("  [WARN: %llu huge errors >= 0.5]",
			(unsigned long long)h->huge);
	}
	printf("\n");
	/* RMSE over all checkpoints, including huge and tiny ones. */
	if (calls != 0 && h->sum_sq > 0) {
		double rmse = sqrt(h->sum_sq / (double)calls);
		printf("  rmse %s: %.6e (e=%.2f) over %llu calls\n",
			name, rmse, -log2(rmse), (unsigned long long)calls);
	}
	if (n != 0) {
		printf("  hist %s:", name);
		for (int e = 1; e <= PREC_EVAL_E_MAX; e ++) {
			if (h->hist[e] != 0) {
				printf(" %d:%llu", e,
					(unsigned long long)h->hist[e]);
			}
		}
		printf("\n");
	}
}

/* Per-key RMSE exponents (e = -log2 rmse) and worst calls. */
typedef struct {
	uint64_t keys;
	double e_min, e_p001, e_p01, e_median, e_max;
	uint64_t worst_key;
	int w_min, w_median;
	uint64_t w_key;
} perkey_stats;

static int
cmp_double(const void *a, const void *b)
{
	double x = *(const double *)a, y = *(const double *)b;
	return (x > y) - (x < y);
}

static int
cmp_int(const void *a, const void *b)
{
	int x = *(const int *)a, y = *(const int *)b;
	return (x > y) - (x < y);
}

static void
perkey_reduce(const perkey_rec *recs, size_t n, int want_sigma,
	perkey_stats *ps)
{
	double *ev = malloc(n * sizeof *ev);
	int *wv = malloc(n * sizeof *wv);

	memset(ps, 0, sizeof *ps);
	ps->keys = n;
	if (n == 0 || ev == 0 || wv == 0) {
		free(ev);
		free(wv);
		return;
	}
	ps->e_min = PREC_EVAL_E_MAX + 1;
	ps->w_min = PREC_EVAL_E_MAX + 1;
	for (size_t i = 0; i < n; i ++) {
		double sq = want_sigma ? recs[i].sq_sig : recs[i].sq_mu;
		int w = want_sigma ? recs[i].min_sig : recs[i].min_mu;
		double e = PREC_EVAL_E_MAX + 1;
		if (recs[i].calls != 0 && sq > 0) {
			e = -log2(sqrt(sq / (double)recs[i].calls));
		}
		ev[i] = e;
		wv[i] = w;
		if (e < ps->e_min) {
			ps->e_min = e;
			ps->worst_key = recs[i].key;
		}
		if (w < ps->w_min) {
			ps->w_min = w;
			ps->w_key = recs[i].key;
		}
	}
	qsort(ev, n, sizeof *ev, cmp_double);
	qsort(wv, n, sizeof *wv, cmp_int);
	ps->e_p001 = ev[n / 1000];
	ps->e_p01 = ev[n / 100];
	ps->e_median = ev[n / 2];
	ps->e_max = ev[n - 1];
	ps->w_median = wv[n / 2];
	free(ev);
	free(wv);
}

static void
print_perkey_summary(const char *name, const perkey_stats *ps)
{
	if (ps->keys == 0) {
		return;
	}
	printf("perkey %-10s keys=%llu  rmse-e: min(worst)=%.3f (key %llu)"
		" p0.1=%.3f p1=%.3f median=%.3f max=%.3f"
		"  worst-call: min=%d (key %llu) median=%d\n",
		name, (unsigned long long)ps->keys, ps->e_min,
		(unsigned long long)ps->worst_key, ps->e_p001, ps->e_p01,
		ps->e_median, ps->e_max, ps->w_min,
		(unsigned long long)ps->w_key, ps->w_median);
}

/* Everything a summary / sidecar carries besides the accumulators. */
typedef struct {
	char target[64], variant[160], master[16], tree[64];
	int metric, frac_bits, ref_bits;
	unsigned logn;
	uint64_t seed, jobs, keys, sigs, keyfirst, signatures;
	uint64_t verify_ok, verify_fail;
} run_meta;

static void
print_summary(const run_meta *m, const prec_eval_stats *st,
	const perkey_stats *pk_mu, const perkey_stats *pk_sig)
{
	uint64_t attempts = st->calls >> (m->logn + 1);

	/* mu is absolute and sigma relative for every target; 'class'
	   only says fixed-point or float. */
	printf("== precision evaluation: target=%s variant=%s logn=%u"
		" runs=%llu ref_bits=%d master=%s class=%s",
		m->target, m->variant, m->logn,
		(unsigned long long)m->signatures, m->ref_bits,
		m->master, (m->metric == PREC_TGT_METRIC_FXP) ? "fxp" : "float");
	if (m->metric == PREC_TGT_METRIC_FXP) {
		printf(" frac_bits=%d", m->frac_bits);
	}
	printf(" ==\n");
	printf("setup: keys=%llu sigs=%llu keyfirst=%llu seed=%llu"
		" jobs=%llu tree=%s\n",
		(unsigned long long)m->keys, (unsigned long long)m->sigs,
		(unsigned long long)m->keyfirst, (unsigned long long)m->seed,
		(unsigned long long)m->jobs, m->tree);
	printf("verify: ok=%llu fail=%llu   sampler calls=%llu"
		"   attempts=%llu (retries=%llu)\n",
		(unsigned long long)m->verify_ok,
		(unsigned long long)m->verify_fail,
		(unsigned long long)st->calls,
		(unsigned long long)attempts,
		(unsigned long long)(attempts - m->signatures));
	printf("rounding divergence (ref vs target): "
		"rint %llu/%llu  floor %llu/%llu  trunc %llu/%llu"
		"  p63 %llu/%llu\n",
		(unsigned long long)st->rint_div,
		(unsigned long long)st->rint_calls,
		(unsigned long long)st->floor_div,
		(unsigned long long)st->floor_calls,
		(unsigned long long)st->trunc_div,
		(unsigned long long)st->trunc_calls,
		(unsigned long long)st->p63_div,
		(unsigned long long)st->p63_calls);
	print_hist_summary("mu_abs", &st->mu, st->calls);
	print_hist_summary(sigma_label(), &st->sigma, st->calls);
	if (pk_mu != 0 && pk_sig != 0) {
		print_perkey_summary("mu_abs", pk_mu);
		print_perkey_summary(sigma_label(), pk_sig);
	}
}

/* ==================================================================== */
/* Stats sidecar (additive, mergeable).                                 */

static void
write_hist(FILE *f, const char *name, const prec_eval_hist *h)
{
	fprintf(f, "hist %s %llu %llu %llu", name,
		(unsigned long long)h->exact, (unsigned long long)h->huge,
		(unsigned long long)h->tiny);
	for (int e = 1; e <= PREC_EVAL_E_MAX; e ++) {
		if (h->hist[e] != 0) {
			fprintf(f, " %d:%llu", e,
				(unsigned long long)h->hist[e]);
		}
	}
	fprintf(f, "\n");
}

static void
write_perkey_stats(FILE *f, const char *name, const perkey_stats *ps)
{
	if (ps->keys == 0) {
		return;
	}
	fprintf(f, "perkey %s %llu %.6f %llu %.6f %.6f %.6f %.6f %d %llu %d\n",
		name, (unsigned long long)ps->keys, ps->e_min,
		(unsigned long long)ps->worst_key, ps->e_p001, ps->e_p01,
		ps->e_median, ps->e_max, ps->w_min,
		(unsigned long long)ps->w_key, ps->w_median);
}

static int
write_stats(const char *path, const run_meta *m, const prec_eval_stats *st,
	const perkey_stats *pk_mu, const perkey_stats *pk_sig)
{
	FILE *f = fopen(path, "w");

	if (f == 0) {
		perror(path);
		return 0;
	}

	fprintf(f, "prec-stats v3\n");
	fprintf(f, "target %s\n", m->target);
	fprintf(f, "variant %s\n", m->variant);
	fprintf(f, "metric %d\n", m->metric);
	fprintf(f, "frac_bits %d\n", m->frac_bits);
	fprintf(f, "ref_bits %d\n", m->ref_bits);
	fprintf(f, "master %s\n", m->master);
	fprintf(f, "tree %s\n", m->tree);
	fprintf(f, "logn %u\n", m->logn);
	fprintf(f, "seed %llu\n", (unsigned long long)m->seed);
	fprintf(f, "jobs %llu\n", (unsigned long long)m->jobs);
	fprintf(f, "keys %llu\n", (unsigned long long)m->keys);
	fprintf(f, "sigs %llu\n", (unsigned long long)m->sigs);
	fprintf(f, "keyfirst %llu\n", (unsigned long long)m->keyfirst);
	/* 'runs' counts signatures. */
	fprintf(f, "runs %llu\n", (unsigned long long)m->signatures);
	fprintf(f, "verify %llu %llu\n", (unsigned long long)m->verify_ok,
		(unsigned long long)m->verify_fail);
	fprintf(f, "calls %llu\n", (unsigned long long)st->calls);
	fprintf(f, "rint %llu %llu\n", (unsigned long long)st->rint_calls,
		(unsigned long long)st->rint_div);
	fprintf(f, "floor %llu %llu\n", (unsigned long long)st->floor_calls,
		(unsigned long long)st->floor_div);
	fprintf(f, "trunc %llu %llu\n", (unsigned long long)st->trunc_calls,
		(unsigned long long)st->trunc_div);
	fprintf(f, "p63 %llu %llu\n", (unsigned long long)st->p63_calls,
		(unsigned long long)st->p63_div);
	write_hist(f, "mu", &st->mu);
	write_hist(f, "sigma", &st->sigma);
	/* %a keeps the double exact through the text file. */
	fprintf(f, "sumsq mu %a\n", st->mu.sum_sq);
	fprintf(f, "sumsq sigma %a\n", st->sigma.sum_sq);
	if (pk_mu != 0) {
		write_perkey_stats(f, "mu", pk_mu);
	}
	if (pk_sig != 0) {
		write_perkey_stats(f, "sigma", pk_sig);
	}
	fclose(f);
	return 1;
}

/* ==================================================================== */
/* Merge mode.                                                          */

static void
merge_hist(prec_eval_hist *h, char *line)
{
	unsigned long long a, b, c;
	char *p = line;

	if (sscanf(p, "%llu %llu %llu", &a, &b, &c) != 3) {
		return;
	}
	h->exact += a;
	h->huge += b;
	h->tiny += c;
	for (p = strchr(p, ':'); p != 0; p = strchr(p + 1, ':')) {
		char *q = p;
		int e;
		unsigned long long n;
		while (q > line && q[-1] != ' ') {
			q --;
		}
		if (sscanf(q, "%d:%llu", &e, &n) == 2
			&& e >= 1 && e <= PREC_EVAL_E_MAX)
		{
			h->hist[e] += n;
		}
	}
}

/* Take a metadata value from a file, refusing one that conflicts with
   an earlier file. "?" matches anything. */
static int
meta_adopt(const char *path, const char *key, char *have, size_t have_len,
	const char *seen)
{
	if (strcmp(have, "?") != 0 && strcmp(have, seen) != 0) {
		fprintf(stderr, "%s: %s '%s' conflicts with '%s'\n",
			path, key, seen, have);
		return 0;
	}
	snprintf(have, have_len, "%s", seen);
	return 1;
}

static int
meta_adopt_u64(const char *path, const char *key, uint64_t *have,
	int *seen_any, uint64_t seen)
{
	if (*seen_any && *have != seen) {
		fprintf(stderr, "%s: %s %llu conflicts with %llu\n",
			path, key, (unsigned long long)seen,
			(unsigned long long)*have);
		return 0;
	}
	*have = seen;
	*seen_any = 1;
	return 1;
}

static int
cmp_perkey(const void *a, const void *b)
{
	uint64_t x = ((const perkey_rec *)a)->key;
	uint64_t y = ((const perkey_rec *)b)->key;
	return (x > y) - (x < y);
}

/* Read the per-key partial records next to a stats file (same path
   with .perkey), appending them to a growable array. */
static int
read_perkey(const char *stats_path, perkey_rec **recs, size_t *n,
	size_t *cap)
{
	char path[1024];
	size_t len = strlen(stats_path);
	FILE *f;
	char line[256];

	if (len < 6 || strcmp(stats_path + len - 6, ".stats") != 0) {
		return 1;
	}
	snprintf(path, sizeof path, "%.*s.perkey", (int)(len - 6),
		stats_path);
	f = fopen(path, "r");
	if (f == 0) {
		return 1;
	}
	while (fgets(line, sizeof line, f) != 0) {
		perkey_rec r;
		unsigned long long k, s, c;
		if (sscanf(line, "%llu %llu %llu %la %la %d %d",
			&k, &s, &c, &r.sq_mu, &r.sq_sig,
			&r.min_mu, &r.min_sig) != 7)
		{
			continue;
		}
		r.key = k;
		r.sigs = s;
		r.calls = c;
		if (*n == *cap) {
			size_t nc = (*cap == 0) ? 1024 : (*cap * 2);
			perkey_rec *nr = realloc(*recs, nc * sizeof *nr);
			if (nr == 0) {
				fclose(f);
				return 0;
			}
			*recs = nr;
			*cap = nc;
		}
		(*recs)[(*n) ++] = r;
	}
	fclose(f);
	return 1;
}

/* Fold records with the same key (a key split across two shards). */
static size_t
fold_perkey(perkey_rec *recs, size_t n)
{
	size_t w = 0;

	if (n == 0) {
		return 0;
	}
	qsort(recs, n, sizeof *recs, cmp_perkey);
	for (size_t i = 1; i < n; i ++) {
		if (recs[i].key == recs[w].key) {
			recs[w].sigs += recs[i].sigs;
			recs[w].calls += recs[i].calls;
			recs[w].sq_mu += recs[i].sq_mu;
			recs[w].sq_sig += recs[i].sq_sig;
			if (recs[i].min_mu < recs[w].min_mu) {
				recs[w].min_mu = recs[i].min_mu;
			}
			if (recs[i].min_sig < recs[w].min_sig) {
				recs[w].min_sig = recs[i].min_sig;
			}
		} else {
			recs[++ w] = recs[i];
		}
	}
	return w + 1;
}

static int
do_merge(int nfiles, char *files[], const char *out)
{
	static prec_eval_stats st;
	run_meta m;
	char sbuf[160];
	int seen_seed = 0, seen_sigs = 0, seen_keyfirst = 0;
	uint64_t keys_sum = 0;
	perkey_rec *recs = 0;
	size_t nrec = 0, cap = 0;
	perkey_stats pk_mu, pk_sig;

	memset(&m, 0, sizeof m);
	snprintf(m.target, sizeof m.target, "?");
	snprintf(m.variant, sizeof m.variant, "?");
	snprintf(m.master, sizeof m.master, "?");
	snprintf(m.tree, sizeof m.tree, "?");
	m.sigs = 1;
	for (int i = 0; i < nfiles; i ++) {
		FILE *f = fopen(files[i], "r");
		char line[16384];

		if (f == 0) {
			perror(files[i]);
			return 1;
		}
		if (fgets(line, sizeof line, f) == 0
			|| (strcmp(line, "prec-stats v2\n") != 0
				&& strcmp(line, "prec-stats v3\n") != 0))
		{
			/* v1 files have no squared-error sums. */
			fprintf(stderr, "%s: not a prec-stats v2/v3 file\n",
				files[i]);
			fclose(f);
			return 1;
		}
		m.jobs ++;
		while (fgets(line, sizeof line, f) != 0) {
			unsigned long long a, b;
			unsigned u;
			double d;
			int ok = 1;
			if (sscanf(line, "metric %d", &m.metric) == 1
				|| sscanf(line, "frac_bits %d",
					&m.frac_bits) == 1
				|| sscanf(line, "ref_bits %d",
					&m.ref_bits) == 1)
			{
				continue;
			}
			if (sscanf(line, "target %63s", sbuf) == 1) {
				ok = meta_adopt(files[i], "target",
					m.target, sizeof m.target, sbuf);
			} else if (sscanf(line, "variant %159s", sbuf) == 1) {
				ok = meta_adopt(files[i], "variant",
					m.variant, sizeof m.variant, sbuf);
			} else if (sscanf(line, "master %15s", sbuf) == 1) {
				ok = meta_adopt(files[i], "master",
					m.master, sizeof m.master, sbuf);
			} else if (sscanf(line, "tree %63s", sbuf) == 1) {
				ok = meta_adopt(files[i], "tree",
					m.tree, sizeof m.tree, sbuf);
			} else if (sscanf(line, "logn %u", &u) == 1) {
				if (m.logn != 0 && m.logn != u) {
					fprintf(stderr, "%s: logn %u"
						" conflicts with %u\n",
						files[i], u, m.logn);
					ok = 0;
				}
				m.logn = u;
			} else if (sscanf(line, "seed %llu", &a) == 1) {
				ok = meta_adopt_u64(files[i], "seed",
					&m.seed, &seen_seed, a);
			} else if (sscanf(line, "sigs %llu", &a) == 1) {
				ok = meta_adopt_u64(files[i], "sigs",
					&m.sigs, &seen_sigs, a);
			} else if (sscanf(line, "keyfirst %llu", &a) == 1) {
				ok = meta_adopt_u64(files[i], "keyfirst",
					&m.keyfirst, &seen_keyfirst, a);
			} else if (sscanf(line, "keys %llu", &a) == 1) {
				keys_sum += a;
			} else if (sscanf(line, "runs %llu", &a) == 1) {
				m.signatures += a;
			} else if (sscanf(line, "verify %llu %llu",
				&a, &b) == 2)
			{
				m.verify_ok += a;
				m.verify_fail += b;
			} else if (sscanf(line, "calls %llu", &a) == 1) {
				st.calls += a;
			} else if (sscanf(line, "rint %llu %llu",
				&a, &b) == 2)
			{
				st.rint_calls += a;
				st.rint_div += b;
			} else if (sscanf(line, "floor %llu %llu",
				&a, &b) == 2)
			{
				st.floor_calls += a;
				st.floor_div += b;
			} else if (sscanf(line, "trunc %llu %llu",
				&a, &b) == 2)
			{
				st.trunc_calls += a;
				st.trunc_div += b;
			} else if (sscanf(line, "p63 %llu %llu",
				&a, &b) == 2)
			{
				st.p63_calls += a;
				st.p63_div += b;
			} else if (strncmp(line, "hist mu ", 8) == 0) {
				merge_hist(&st.mu, line + 8);
			} else if (strncmp(line, "hist sigma ", 11) == 0) {
				merge_hist(&st.sigma, line + 11);
			} else if (sscanf(line, "sumsq mu %la", &d) == 1) {
				st.mu.sum_sq += d;
			} else if (sscanf(line, "sumsq sigma %la", &d) == 1) {
				st.sigma.sum_sq += d;
			}
			if (!ok) {
				fclose(f);
				return 1;
			}
		}
		fclose(f);
		if (!read_perkey(files[i], &recs, &nrec, &cap)) {
			fprintf(stderr, "%s: out of memory on perkey\n",
				files[i]);
			return 1;
		}
	}
	nrec = fold_perkey(recs, nrec);
	/* Distinct keys from the fold if available. */
	m.keys = (nrec != 0) ? nrec : keys_sum;
	perkey_reduce(recs, nrec, 0, &pk_mu);
	perkey_reduce(recs, nrec, 1, &pk_sig);
	if (out != 0) {
		char path[1024];
		snprintf(path, sizeof path, "%s.stats", out);
		if (!write_stats(path, &m, &st, &pk_mu, &pk_sig)) {
			return 1;
		}
		if (nrec != 0) {
			FILE *f;
			snprintf(path, sizeof path, "%s.perkey.csv", out);
			f = fopen(path, "w");
			if (f == 0) {
				perror(path);
				return 1;
			}
			fprintf(f, "key,sigs,calls,rmse_mu_abs,rmse_%s,"
				"e_mu_abs,e_%s,min_mu_abs_e,min_%s_e\n",
				sigma_label(), sigma_label(), sigma_label());
			for (size_t i = 0; i < nrec; i ++) {
				const perkey_rec *r = &recs[i];
				double rm = 0, rs = 0;
				if (r->calls != 0) {
					rm = sqrt(r->sq_mu / (double)r->calls);
					rs = sqrt(r->sq_sig / (double)r->calls);
				}
				fprintf(f, "%llu,%llu,%llu,%.9e,%.9e,"
					"%.3f,%.3f,%d,%d\n",
					(unsigned long long)r->key,
					(unsigned long long)r->sigs,
					(unsigned long long)r->calls, rm, rs,
					rm > 0 ? -log2(rm) : 0.0,
					rs > 0 ? -log2(rs) : 0.0,
					r->min_mu, r->min_sig);
			}
			fclose(f);
		}
	}
	free(recs);
	print_summary(&m, &st, &pk_mu, &pk_sig);
	return 0;
}

/* ==================================================================== */

int
main(int argc, char *argv[])
{
	unsigned logn = 10;
	uint64_t keys = 25, sigs = 1, keyfirst = 0, seed = 1;
	unsigned shard = 0, jobs = 1;
	int want_csv = 0;
	const char *out = 0;

	for (int i = 1; i < argc; i ++) {
		if (strcmp(argv[i], "--merge") == 0) {
			return do_merge(argc - i - 1, argv + i + 1, out);
		} else if (strcmp(argv[i], "--variant-tokens") == 0) {
			printf("%s\n", name_tokens());
			return 0;
		} else if (strcmp(argv[i], "--csv-header") == 0) {
			printf("run,mu_abs_e,%s_e\n", sigma_label());
			return 0;
		} else if (strcmp(argv[i], "--worst-header") == 0) {
			printf("key,sig,min_mu_abs_e,min_%s_e\n",
				sigma_label());
			return 0;
		} else if (strcmp(argv[i], "--csv") == 0) {
			want_csv = 1;
		} else if (strcmp(argv[i], "--logn") == 0 && i + 1 < argc) {
			logn = (unsigned)strtoul(argv[++ i], 0, 10);
		} else if (strcmp(argv[i], "--runs") == 0 && i + 1 < argc) {
			keys = strtoull(argv[++ i], 0, 10);
		} else if (strcmp(argv[i], "--sigs") == 0 && i + 1 < argc) {
			sigs = strtoull(argv[++ i], 0, 10);
		} else if (strcmp(argv[i], "--key-first") == 0
			&& i + 1 < argc)
		{
			keyfirst = strtoull(argv[++ i], 0, 10);
		} else if (strcmp(argv[i], "--seed") == 0 && i + 1 < argc) {
			seed = strtoull(argv[++ i], 0, 10);
		} else if (strcmp(argv[i], "--tree") == 0 && i + 1 < argc) {
			tree_id = argv[++ i];
		} else if (strcmp(argv[i], "--shard") == 0 && i + 1 < argc) {
			shard = (unsigned)strtoul(argv[++ i], 0, 10);
		} else if (strcmp(argv[i], "--jobs") == 0 && i + 1 < argc) {
			jobs = (unsigned)strtoul(argv[++ i], 0, 10);
		} else if (strcmp(argv[i], "--out") == 0 && i + 1 < argc) {
			out = argv[++ i];
		} else {
			fprintf(stderr, "usage: eval_precision"
				" [--logn L] [--runs KEYS] [--sigs M]"
				" [--key-first F] [--seed S] [--tree ID]"
				" [--shard I --jobs J] [--out PREFIX]"
				" [--csv]\n"
				"       eval_precision [--out PREFIX]"
				" [--tree ID] --merge FILE.stats...\n"
				"       eval_precision --csv-header\n"
				"       eval_precision --worst-header\n"
				"       eval_precision --variant-tokens\n");
			return 2;
		}
	}
	if (logn < 2 || logn > 10 || jobs == 0 || shard >= jobs
		|| keys == 0 || sigs == 0)
	{
		fprintf(stderr, "eval_precision: bad parameters\n");
		return 2;
	}

	fndsa_prec_init();
	fndsa_prec_stats_reset();

	FILE *worst = 0, *perkey = 0;
	char path[1024];
	if (out != 0) {
		snprintf(path, sizeof path, "%s.worst", out);
		worst = fopen(path, "w");
		if (worst == 0) {
			perror(path);
			return 1;
		}
		snprintf(path, sizeof path, "%s.perkey", out);
		perkey = fopen(path, "w");
		if (perkey == 0) {
			perror(path);
			return 1;
		}
		if (want_csv) {
			snprintf(path, sizeof path, "%s.csv", out);
			csv = fopen(path, "w");
			if (csv == 0) {
				perror(path);
				return 1;
			}
		}
	}
	fndsa_prec_set_sink(eval_sink, 0);
	signal(SIGALRM, alarm_handler);

	/* Signature index g in [first, last), split evenly across shards:
	   key r = keyfirst + g / sigs, signature v = g % sigs. Seeds depend
	   only on (seed, logn, r, v), so results do not depend on the
	   number of shards. */
	uint64_t total = keys * sigs;
	uint64_t first = total * shard / jobs;
	uint64_t last = total * (shard + 1) / jobs;
	uint64_t verify_ok = 0, verify_fail = 0, shard_keys = 0;
	uint64_t cur_key = 0;
	int have_key = 0;
	perkey_rec pk;
	const prec_eval_stats *pst = fndsa_prec_stats();

	uint8_t *sk = malloc(FNDSA_SIGN_KEY_SIZE(10));
	uint8_t *vk = malloc(FNDSA_VRFY_KEY_SIZE(10));
	uint8_t *sig = malloc(FNDSA_SIGNATURE_SIZE(10));

	memset(&pk, 0, sizeof pk);
	for (uint64_t g = first; g < last; g ++) {
		uint64_t r = keyfirst + g / sigs, v = g % sigs;
		char kseed[96], sseed[96], msg[96];
		size_t siglen;
		uint64_t c0 = pst->calls;
		double m0 = pst->mu.sum_sq, s0 = pst->sigma.sum_sq;

		if (!have_key || r != cur_key) {
			if (have_key && perkey != 0) {
				perkey_write(perkey, &pk);
			}
			snprintf(kseed, sizeof kseed,
				"prec-keygen-%llu-%u-%llu",
				(unsigned long long)seed, logn,
				(unsigned long long)r);
			alarm(600);
			fndsa_keygen_seeded(logn, kseed, strlen(kseed), sk, vk);
			alarm(0);
			memset(&pk, 0, sizeof pk);
			pk.key = r;
			pk.min_mu = PREC_EVAL_E_MAX + 1;
			pk.min_sig = PREC_EVAL_E_MAX + 1;
			cur_key = r;
			have_key = 1;
			shard_keys ++;
		}
		if (sigs == 1) {
			snprintf(sseed, sizeof sseed,
				"prec-sign-%llu-%u-%llu",
				(unsigned long long)seed, logn,
				(unsigned long long)r);
			snprintf(msg, sizeof msg,
				"prec-message-%llu-%u-%llu",
				(unsigned long long)seed, logn,
				(unsigned long long)r);
		} else {
			snprintf(sseed, sizeof sseed,
				"prec-sign-%llu-%u-%llu-%llu",
				(unsigned long long)seed, logn,
				(unsigned long long)r, (unsigned long long)v);
			snprintf(msg, sizeof msg,
				"prec-message-%llu-%u-%llu-%llu",
				(unsigned long long)seed, logn,
				(unsigned long long)r, (unsigned long long)v);
		}
		cur_run = g;
		run_worst_reset();
		alarm(600);
		siglen = fndsa_sign_seeded(sk, FNDSA_SIGN_KEY_SIZE(logn),
			NULL, 0, FNDSA_HASH_ID_RAW, msg, strlen(msg),
			sseed, strlen(sseed),
			sig, FNDSA_SIGNATURE_SIZE(logn));
		alarm(0);
		if (siglen != 0 && fndsa_verify(sig, siglen,
			vk, FNDSA_VRFY_KEY_SIZE(logn),
			NULL, 0, FNDSA_HASH_ID_RAW, msg, strlen(msg)))
		{
			verify_ok ++;
		} else {
			verify_fail ++;
		}
		pk.sigs ++;
		pk.calls += pst->calls - c0;
		pk.sq_mu += pst->mu.sum_sq - m0;
		pk.sq_sig += pst->sigma.sum_sq - s0;
		if (run_min_mu < pk.min_mu) {
			pk.min_mu = run_min_mu;
		}
		if (run_min_sig < pk.min_sig) {
			pk.min_sig = run_min_sig;
		}
		if (worst != 0) {
			fprintf(worst, "%llu,%llu,%d,%d\n",
				(unsigned long long)r, (unsigned long long)v,
				run_min_mu, run_min_sig);
		}
	}
	if (have_key && perkey != 0) {
		perkey_write(perkey, &pk);
	}
	free(sk);
	free(vk);
	free(sig);

	if (csv != 0) {
		fclose(csv);
		csv = 0;
	}
	if (worst != 0) {
		fclose(worst);
	}
	if (perkey != 0) {
		fclose(perkey);
	}
	run_meta m;
	memset(&m, 0, sizeof m);
	snprintf(m.target, sizeof m.target, "%s", fndsa_prec_target());
	snprintf(m.variant, sizeof m.variant, "%s",
		variant_label(name_tokens()));
	snprintf(m.master, sizeof m.master, "%s",
		fndsa_prec_master_is_ref() ? "ref" : "target");
	snprintf(m.tree, sizeof m.tree, "%s", tree_id);
	m.metric = fndsa_prec_metric();
	m.frac_bits = fndsa_prec_frac_bits();
	m.ref_bits = fndsa_prec_ref_bits();
	m.logn = logn;
	m.seed = seed;
	m.jobs = 1;
	m.keys = shard_keys;
	m.sigs = sigs;
	m.keyfirst = keyfirst;
	m.signatures = last - first;
	m.verify_ok = verify_ok;
	m.verify_fail = verify_fail;
	if (out != 0) {
		snprintf(path, sizeof path, "%s.stats", out);
		if (!write_stats(path, &m, pst, 0, 0)) {
			return 1;
		}
	}
	if (jobs == 1) {
		print_summary(&m, pst, 0, 0);
	}
	return verify_fail == 0 ? 0 : 1;
}
