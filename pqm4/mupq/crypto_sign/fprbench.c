/*
 * Per-operation cycle benchmark for the fpr backend, a local addition to
 * pqm4. It times the fpr contract ops and the fpoly/FFT routines, so the
 * rows match across backends and implementations.
 *
 * BENCH (fpr_bench_cm4.S) reads the DWT cycle counter around one call.
 * Each figure is the minimum over REPS runs minus a calibration: single
 * fpr ops go through wrappers and subtract a baseline wrapper of the same
 * shape (fpr_bench_wrap.c); fpoly routines subtract the bare harness.
 * Built by 'make pqm4-bench'.
 */

#include "api.h"
#include "hal.h"
#include "sendfn.h"

#include <stdio.h>
#include <stdint.h>
#include <string.h>

/* The fpr contract and fpoly_* routines of this directory's backend. */
#include "sign_inner.h"

/* ---------------------------------------------------------------------- */
/* Benchmark entry points (fpr_bench_cm4.S): each returns the cycle count of
   one call to its target. */

extern uint32_t bench_none(void);
extern uint32_t bench_base0(void);
extern uint32_t bench_base1(void);
extern uint32_t bench_base1nc(void);
extern uint32_t bench_base2(void);
extern uint32_t bench_base2nc(void);
extern uint32_t bench_base2s2(void);
extern uint32_t bench_base2s2nc(void);
extern uint32_t bench_base_addsub(void);
extern uint32_t bench_base_recip_x(void);
extern uint32_t bench_base_mul_recip(void);
extern uint32_t bench_base_sqrt_sigma(void);

extern uint32_t bench_fpr_add(void);
extern uint32_t bench_fpr_mul(void);
extern uint32_t bench_fpr_div(void);
extern uint32_t bench_fpr_sqrt(void);
extern uint32_t bench_fpr_sqrt_sigma(void);
extern uint32_t bench_fpr_add_sub(void);
extern uint32_t bench_fpr_sub(void);
extern uint32_t bench_fpr_neg(void);
extern uint32_t bench_fpr_half(void);
extern uint32_t bench_fpr_double(void);
extern uint32_t bench_fpr_sqr(void);
extern uint32_t bench_fpr_inv(void);
extern uint32_t bench_fpr_recip_x(void);
extern uint32_t bench_fpr_mul_recip(void);
extern uint32_t bench_fpr_div2e(void);
extern uint32_t bench_fpr_mul2e(void);
extern uint32_t bench_fpr_of(void);
extern uint32_t bench_fpr_of32(void);
extern uint32_t bench_fpr_rint(void);
extern uint32_t bench_fpr_floor(void);
extern uint32_t bench_fpr_trunc(void);
extern uint32_t bench_fpc_mul(void);
extern uint32_t bench_fpc_mul_gm(void);
extern uint32_t bench_fpc_mul_gmj(void);
extern uint32_t bench_fpr_mul_gm(void);

extern uint32_t bench_fpoly_FFT(unsigned logn, fpr *f);
#if FNDSA_FFT_NEW
extern uint32_t bench_fpoly_FFT_int32(unsigned logn, fpr *d,
	const int32_t *a);
extern uint32_t bench_fpoly_FFT_int32_selfadj(unsigned logn, fpr *d,
	const int32_t *a);
extern uint32_t bench_fpoly_FFT_int8(unsigned logn, fpr *d,
	const int8_t *a);
extern uint32_t bench_fpoly_FFT_uint16(unsigned logn, fpr *d,
	const uint16_t *a);
extern uint32_t bench_expand_int8(unsigned logn, fpr *d, const int8_t *a);
extern uint32_t bench_expand_uint16(unsigned logn, fpr *d,
	const uint16_t *a);
extern uint32_t bench_expand_uint16_fpr(unsigned logn, fpr *d,
	const uint16_t *a);
#endif
extern uint32_t bench_fpoly_iFFT(unsigned logn, fpr *f);
extern uint32_t bench_fpoly_set_small(unsigned logn, fpr *d, const int8_t *f);
extern uint32_t bench_fpoly_add(unsigned logn, fpr *a, const fpr *b);
extern uint32_t bench_fpoly_sub(unsigned logn, fpr *a, const fpr *b);
extern uint32_t bench_fpoly_neg(unsigned logn, fpr *a);
extern uint32_t bench_fpoly_mul_fft(unsigned logn, fpr *a, const fpr *b);
extern uint32_t bench_fpoly_LDL_fft(unsigned logn,
	const fpr *g00, fpr *g01, fpr *g11);
extern uint32_t bench_fpoly_split_fft(unsigned logn,
	fpr *f0, fpr *f1, const fpr *f);
extern uint32_t bench_fpoly_split_selfadj(unsigned logn,
	fpr *f0, fpr *f1, const fpr *f);
extern uint32_t bench_fpoly_merge_fft(unsigned logn,
	fpr *f, const fpr *f0, const fpr *f1);
extern uint32_t bench_fpoly_gram_fft(unsigned logn);
extern uint32_t bench_fpoly_mulconst(unsigned logn, fpr *a);

/* Operand globals consumed by the wrappers (fpr_bench_wrap.c). */
extern volatile fpr fndsa_bw_x, fndsa_bw_y;
extern volatile fpr fndsa_bw_sink_a, fndsa_bw_sink_b;
#define SYSCNT_ADDR   0xE0001004u   /* DWT_CYCCNT, as in fpr_bench_cm4.S */
extern volatile int32_t fndsa_bw_i;
extern volatile unsigned fndsa_bw_e;
extern fpr *fndsa_bw_p00, *fndsa_bw_p01, *fndsa_bw_p10, *fndsa_bw_p11;

/* Direct calls to fndsa_fpr_* symbols with operands in registers, as in
   upstream bench_cm4. */
extern const uint32_t fndsa_bench_dir_present[6];
extern uint32_t bench_direct_mul(fpr x, fpr y);
extern uint32_t bench_direct_sqr(fpr x);
extern uint32_t bench_direct_div(fpr x, fpr y);
extern uint32_t bench_direct_sqrt(fpr x);
extern uint32_t bench_direct_of32(int32_t i);
extern uint32_t bench_direct_scaled(int64_t i, int32_t sc);

/* ---------------------------------------------------------------------- */

static void
enable_cyccnt(void)
{
	volatile uint32_t *DEMCR      = (volatile uint32_t *)0xE000EDFC;
	volatile uint32_t *DWT_CTRL   = (volatile uint32_t *)0xE0001000;
	volatile uint32_t *DWT_CYCCNT = (volatile uint32_t *)0xE0001004;
	volatile uint32_t *DWT_LAR    = (volatile uint32_t *)0xE0001FB0;

	*DEMCR |= 0x01000000;    /* TRCENA */
	*DWT_LAR = 0xC5ACCE55;   /* unlock */
	*DWT_CYCCNT = 0;
	*DWT_CTRL |= 1;          /* CYCCNTENA */
}

#define REPS   16

static uint32_t cal_direct;   /* bare harness (BENCH around 'bx lr') */
static uint32_t cal_base0;    /* harness + wrapper with a nop inner call */
static uint32_t cal_base1;    /* + one operand load, inner call, a store */
static uint32_t cal_base1nc;  /* one load + store, NO inner call (inline ops) */
static uint32_t cal_base2;    /* + two operand loads, inner call, a store */
static uint32_t cal_base2nc;  /* two loads + one store, NO inner call */
static uint32_t cal_base2s2;  /* two loads, inner call, two stores */
static uint32_t cal_base2s2nc; /* two loads + two stores, NO inner call */
static uint32_t cal_addsub;   /* the FPR_ADD_SUB call-site shape */
static uint32_t cal_recip_x;  /* the FPR_RECIP_X expansion shape */
static uint32_t cal_mul_recip; /* the FPR_MUL_RECIP expansion shape */
static uint32_t cal_sqrt_sigma; /* the FPR_SQRT_MUL_INV_SIGMA shape */

/* Measurements run with interrupts masked (the HAL's SysTick would land
   in long runs) and after one discarded warm-up run (cold caches). */
#define IRQ_OFF()   __asm__ volatile("cpsid i" ::: "memory")
#define IRQ_ON()    __asm__ volatile("cpsie i" ::: "memory")

/* Spread (max - min) of the last min_of run. */
static uint32_t g_spread;

static uint32_t
min_of(uint32_t (*fn)(void))
{
	uint32_t best = 0xFFFFFFFF, worst = 0;
	IRQ_OFF();
	(void)fn();
	for (int i = 0; i < REPS; i++) {
		uint32_t c = fn();
		if (c < best) {
			best = c;
		}
		if (c > worst) {
			worst = c;
		}
	}
	IRQ_ON();
	g_spread = worst - best;
	return best;
}

static void
report_j(const char *name, uint32_t raw, uint32_t cal, uint32_t spread)
{
	char line[64];
	uint32_t v = (raw > cal) ? (raw - cal) : 0;
	sprintf(line, "%-20s %6lu  spread %lu", name,
		(unsigned long)v, (unsigned long)spread);
	hal_send_str(line);
}

static void
report(const char *name, uint32_t raw, uint32_t cal)
{
	report_j(name, raw, cal, g_spread);
}

/* A scalar benchmark: min over REPS, reported against its baseline. */
#define SCALAR(name, fn, cal)   report((name), min_of(&(fn)), (cal))

/* ---------------------------------------------------------------------- */
/* Working buffers for the largest degree; static to spare the stack. */

#define BENCH_LOGN_MAX   10
#define BENCH_N_MAX      ((size_t)1 << BENCH_LOGN_MAX)

static int8_t small_src[BENCH_N_MAX];
static int32_t int_src[BENCH_N_MAX];
static uint16_t uint16_src[BENCH_N_MAX];
static fpr fft_input[BENCH_N_MAX];
static fpr pristine[BENCH_N_MAX];   /* an FFT-domain polynomial, kept intact */
static fpr work[BENCH_N_MAX];       /* destination for destructive routines */
static fpr operand[BENCH_N_MAX];    /* second operand */
/* Basis B = [[b00,b01],[b10,b11]] in FFT domain and a scratch copy,
   re-seeded before each run since fpoly_gram_fft overwrites it. */
static fpr bb00[BENCH_N_MAX], bb01[BENCH_N_MAX];
static fpr bb10[BENCH_N_MAX], bb11[BENCH_N_MAX];
static fpr w00[BENCH_N_MAX], w01[BENCH_N_MAX], w10[BENCH_N_MAX];

#if FNDSA_FFT_NEW
/* Input widening for the compact FFT rows (benchmark only). */
void
fndsa_bench_expand_int8(unsigned logn, fpr *d, const int8_t *a)
{
	size_t n = (size_t)1 << logn;
	int32_t *cc = (int32_t *)(void *)(d + n) - n;
	for (size_t i = 0; i < n; i ++) {
		cc[i] = a[i];
	}
}

void
fndsa_bench_expand_uint16(unsigned logn, fpr *d, const uint16_t *a)
{
	size_t n = (size_t)1 << logn;
	int32_t *cc = (int32_t *)(void *)(d + n) - n;
	for (size_t i = 0; i < n; i ++) {
		cc[i] = a[i];
	}
}

void
fndsa_bench_expand_uint16_fpr(unsigned logn, fpr *d, const uint16_t *a)
{
	size_t n = (size_t)1 << logn;
	for (size_t i = 0; i < n; i ++) {
		d[i] = fpr_of(a[i]);
	}
}
#endif

static void
fill_small(int seed)
{
	for (size_t i = 0; i < BENCH_N_MAX; i++) {
		int v = (int)((i * (size_t)(3 + 2 * seed) + (size_t)seed) % 11) - 5;
		small_src[i] = (int8_t)(v == 0 ? 1 + seed : v);
	}
}

/* Same input for the classic and the compact FFT. The compact source is
   placed in its workspace before timing, as in basis_to_FFT(). */
static void
setup_fft_input(unsigned logn)
{
	size_t n = (size_t)1 << logn;
	for (size_t i = 0; i < n; i ++) {
		int32_t x = (int32_t)((i * 13 + 17) % 251) - 125;
		int_src[i] = x;
		uint16_src[i] = (uint16_t)((i * 73 + 29) % 12289);
		fft_input[i] = fpr_of32(x);
	}
}

static void
reseed_fft_old(unsigned logn)
{
	size_t n = (size_t)1 << logn;
	memcpy(work, fft_input, n * sizeof *work);
}

#if FNDSA_FFT_NEW
static void
reseed_fft_int32(unsigned logn)
{
	size_t n = (size_t)1 << logn;
	int32_t *cc = (int32_t *)(void *)(work + n) - n;
	memcpy(cc, int_src, n * sizeof *cc);
}

static int32_t *
fft_int32_alias(unsigned logn)
{
	size_t n = (size_t)1 << logn;
	return (int32_t *)(void *)(work + n) - n;
}

static void
reseed_fft_int32_selfadj(unsigned logn)
{
	size_t hn = (size_t)1 << (logn - 1);
	int32_t *aa = (int32_t *)(void *)(work + hn) - hn;
	memcpy(aa, int_src, hn * sizeof *aa);
}

static int32_t *
fft_int32_selfadj_alias(unsigned logn)
{
	size_t hn = (size_t)1 << (logn - 1);
	return (int32_t *)(void *)(work + hn) - hn;
}
#endif

/* Signing-like FFT-domain data, without zeros, so division and square
   root do real work. */
static void
setup(unsigned logn)
{
	size_t n = (size_t)1 << logn;
	size_t sz = n * sizeof(fpr);

	fill_small(1);
	fpoly_set_small(logn, pristine, small_src);
	fpoly_FFT(logn, pristine);

	fill_small(2);
	fpoly_set_small(logn, operand, small_src);
	fpoly_FFT(logn, operand);

	fill_small(3);
	fpoly_set_small(logn, bb00, small_src);
	fpoly_FFT(logn, bb00);
	fill_small(4);
	fpoly_set_small(logn, bb01, small_src);
	fpoly_FFT(logn, bb01);
	fill_small(5);
	fpoly_set_small(logn, bb10, small_src);
	fpoly_FFT(logn, bb10);
	fill_small(6);
	fpoly_set_small(logn, bb11, small_src);
	fpoly_FFT(logn, bb11);

	memcpy(work, pristine, sz);
	memcpy(w00, bb00, sz);
	memcpy(w01, bb01, sz);
	memcpy(w10, bb10, sz);
}

/* Re-seed the buffers a destructive routine overwrites. */
static void
reseed(unsigned logn)
{
	size_t sz = ((size_t)1 << logn) * sizeof(fpr);
	memcpy(work, pristine, sz);
	memcpy(w00, bb00, sz);
	memcpy(w01, bb01, sz);
	memcpy(w10, bb10, sz);
}

/* Gram matrix, as signing builds it: g00 <- w00, g01 <- w01, g11 <- w10. */
static void
reseed_gram(unsigned logn)
{
	reseed(logn);
	fpoly_gram_fft(logn, w00, w01, w10, bb11);
}

/* min over REPS of a poly benchmark, re-seeding between runs. */
#define POLY(name, expr, seedfn)   do { \
		uint32_t best = 0xFFFFFFFF; \
		uint32_t worst = 0; \
		IRQ_OFF(); \
		seedfn(logn); \
		(void)(expr); \
		for (int rep_ = 0; rep_ < REPS; rep_++) { \
			seedfn(logn); \
			uint32_t c_ = (expr); \
			if (c_ < best) { \
				best = c_; \
			} \
			if (c_ > worst) { \
				worst = c_; \
			} \
		} \
		IRQ_ON(); \
		report_j((name), best, cal_direct, worst - best); \
	} while (0)

static void
bench_scalar(void)
{
	hal_send_str("--- fpr operations, one call through a C wrapper ---");

	/* Representative operands: normal range, non-zero, not powers of two. */
	fpr x = fpr_of(1234567);
	fpr y = fpr_of(7654321);
	fndsa_bw_x = fpr_div2e(x, 7);
	fndsa_bw_y = fpr_div2e(y, 3);
	fndsa_bw_i = 12345;
	fndsa_bw_e = 5;

#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264
	/* fxp32.64: only the out-of-line kernels used by signing. The inline
	   ops, sqr and mul_gm are in the direct-call and in-register sections;
	   through this wrapper they would mostly time the wrapper. */
#if FNDSA_FPR_IMPL != FNDSA_FPR_IMPL_CM4
	SCALAR("fpr_mul",     bench_fpr_mul,     cal_base2);
#endif
	SCALAR("FPR_RECIP_X",   bench_fpr_recip_x,   cal_recip_x);
	SCALAR("FPR_MUL_RECIP", bench_fpr_mul_recip, cal_mul_recip);
	SCALAR("FPR_SQRT_SIGMA", bench_fpr_sqrt_sigma, cal_sqrt_sigma);
	SCALAR("FPC_MUL",     bench_fpc_mul,     cal_base2s2);
#if FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4 || FNDSA_TWIDDLE_U64
	SCALAR("FPC_MUL_GM",  bench_fpc_mul_gm,  cal_base2s2);
	SCALAR("FPC_MUL_GMJ", bench_fpc_mul_gmj, cal_base2s2);
#endif
#else
#ifdef FNDSA_FPR_ADD_INLINE
	SCALAR("fpr_add",     bench_fpr_add,     cal_base2nc);
#else
	SCALAR("fpr_add",     bench_fpr_add,     cal_base2);
#endif
#ifdef FNDSA_FPR_SUB_INLINE
	SCALAR("fpr_sub",     bench_fpr_sub,     cal_base2nc);
#else
	SCALAR("fpr_sub",     bench_fpr_sub,     cal_base2);
#endif
#ifdef FNDSA_FPR_ADD_SUB_INLINE
	SCALAR("fpr_add_sub", bench_fpr_add_sub, cal_base2s2nc);
#else
	SCALAR("fpr_add_sub", bench_fpr_add_sub, cal_addsub);
#endif
	SCALAR("fpr_mul",     bench_fpr_mul,     cal_base2);
	SCALAR("fpr_sqr",     bench_fpr_sqr,     cal_base1);
	SCALAR("fpr_div",     bench_fpr_div,     cal_base2);
	SCALAR("fpr_inv",     bench_fpr_inv,     cal_base1);
	/* recip_x runs first and leaves the reciprocal for mul_recip. */
	SCALAR("FPR_RECIP_X",   bench_fpr_recip_x,   cal_recip_x);
	SCALAR("FPR_MUL_RECIP", bench_fpr_mul_recip, cal_mul_recip);
	SCALAR("fpr_sqrt",    bench_fpr_sqrt,    cal_base1);
	SCALAR("FPR_SQRT_SIGMA", bench_fpr_sqrt_sigma, cal_sqrt_sigma);
	SCALAR("FPC_MUL",     bench_fpc_mul,     cal_base2s2);
	SCALAR("fpr_neg",     bench_fpr_neg,     cal_base1nc);
	SCALAR("fpr_half",    bench_fpr_half,    cal_base1nc);
	SCALAR("fpr_double",  bench_fpr_double,  cal_base1nc);
	SCALAR("fpr_div2e",   bench_fpr_div2e,   cal_base1nc);
	SCALAR("fpr_mul2e",   bench_fpr_mul2e,   cal_base1nc);
	SCALAR("fpr_of",      bench_fpr_of,      cal_base1nc);
	SCALAR("fpr_of32",    bench_fpr_of32,    cal_base1);
	SCALAR("fpr_rint",    bench_fpr_rint,    cal_base1nc);
	SCALAR("fpr_floor",   bench_fpr_floor,   cal_base1nc);
	SCALAR("fpr_trunc",   bench_fpr_trunc,   cal_base1nc);
#endif
}

/* Upstream bench_cm4 accounting: one call minus the bench_none
   calibration, as in c-fn-dsa's ttsign.c. Absent targets are skipped. */
static void
bench_direct(void)
{
	uint32_t cal_raw = min_of(&bench_none);
	fpr x = fndsa_bw_x;
	fpr y = fndsa_bw_y;
	uint32_t best, worst, c;

	hal_send_str("--- fpr operations, one direct call from assembly ---");

/* Only fpr_mul has a wrapper for its custom convention; the others need
   sizeof(fpr) <= 8, since a wider fpr is returned through a hidden
   pointer. */
#define DIRECT(idx, name, expr, extra)   do { \
		if (fndsa_bench_dir_present[(idx)] != 0 \
				&& ((idx) == 0 || sizeof(fpr) <= 8)) { \
			best = 0xFFFFFFFF; \
			worst = 0; \
			IRQ_OFF(); \
			(void)(expr); \
			for (int i = 0; i < REPS; i++) { \
				c = (expr); \
				if (c < best) { \
					best = c; \
				} \
				if (c > worst) { \
					worst = c; \
				} \
			} \
			IRQ_ON(); \
			report_j((name), best, cal_raw + (extra), \
				worst - best); \
		} \
	} while (0)

	DIRECT(0, "fpr_mul",     bench_direct_mul(x, y), 0);
	DIRECT(1, "fpr_sqr",     bench_direct_sqr(x), 0);
	DIRECT(2, "fpr_div",     bench_direct_div(x, y), 0);
	DIRECT(3, "fpr_sqrt",    bench_direct_sqrt(x), 0);
	DIRECT(4, "fpr_of32",    bench_direct_of32(12345), 0);
	DIRECT(5, "fpr_scaled",  bench_direct_scaled(12345, -8), 0);
#undef DIRECT
}

/* Inlined ops: INREG_N dependent copies run between two counter reads,
   minus an empty chain. The counter reads take the chain's words as asm
   operands, so the compiler cannot move work across them. Only ops the
   backend inlines are listed; calls are timed in the sections above. */
#define INREG_N   16
#define INREG_W   (sizeof(fpr) / 4)
typedef union { fpr f; uint32_t w[INREG_W]; } inreg_u;

static uint32_t inreg_cal;

#define INREG_TIE(u, t)   do { \
		for (size_t k_ = 0; k_ < INREG_W; k_ ++) { \
			__asm__ volatile("" : "+r" ((u).w[k_]) : "r" (t)); \
		} \
	} while (0)
#define INREG_T0(t, u)   do { \
		__asm__ volatile("ldr %0, [%1]" : "=r" (t) \
			: "r" (SYSCNT_ADDR) : "memory"); \
		INREG_TIE((u), (t)); \
	} while (0)
#define INREG_T1(t, u)   do { \
		uint32_t s_ = 0; \
		for (size_t k_ = 0; k_ < INREG_W; k_ ++) { \
			__asm__ volatile("" : "+r" (s_) : "r" ((u).w[k_])); \
		} \
		__asm__ volatile("ldr %0, [%1]" : "=r" (t) \
			: "r" (SYSCNT_ADDR), "r" (s_) : "memory"); \
	} while (0)
#define INREG_16(...)   do { __VA_ARGS__ __VA_ARGS__ __VA_ARGS__ \
	__VA_ARGS__ __VA_ARGS__ __VA_ARGS__ __VA_ARGS__ __VA_ARGS__ \
	__VA_ARGS__ __VA_ARGS__ __VA_ARGS__ __VA_ARGS__ __VA_ARGS__ \
	__VA_ARGS__ __VA_ARGS__ __VA_ARGS__ } while (0)

/* Makes the compiler treat x as changed, so a chain cannot be cancelled
   or folded (16 negations are the identity). Emits no instructions. */
#define INREG_BAR(u)   do { \
		for (size_t k_ = 0; k_ < INREG_W; k_ ++) { \
			__asm__ volatile("" : "+r" ((u).w[k_])); \
		} \
	} while (0)

/* One measured chain: the statement is applied 16 times to x (and y). */
#define INREG_FN(name, ...) \
	static uint32_t inreg_##name(void) \
	{ \
		inreg_u x, y; \
		uint32_t t0, t1; \
		x.f = fndsa_bw_x; \
		y.f = fndsa_bw_y; \
		INREG_T0(t0, x); \
		INREG_16(__VA_ARGS__); \
		INREG_T1(t1, x); \
		fndsa_bw_sink_a = x.f; \
		fndsa_bw_sink_b = y.f; \
		return t1 - t0; \
	}

INREG_FN(nop,     __asm__ volatile("" : "+r" (x.w[0])); )
#ifdef FNDSA_FPR_ADD_INLINE
INREG_FN(add,     x.f = fpr_add(x.f, y.f); )
#endif
#ifdef FNDSA_FPR_SUB_INLINE
INREG_FN(sub,     x.f = fpr_sub(x.f, y.f); )
#endif
#ifdef FNDSA_FPR_ADD_SUB_INLINE
#define INREG_ADDSUB_STMT   { fpr a_, b_; FPR_ADD_SUB(a_, b_, x.f, y.f); x.f = a_; y.f = b_; }
INREG_FN(add_sub, INREG_ADDSUB_STMT)
#endif
INREG_FN(neg,     x.f = fpr_neg(x.f); INREG_BAR(x); )
INREG_FN(half,    x.f = fpr_half(x.f); )
#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264
INREG_FN(of32,    x.f = fpr_of32((int32_t)(x.w[0] >> 16)); INREG_BAR(x); )
#endif
INREG_FN(rint,    x.w[0] += (uint32_t)fpr_rint(x.f); )
INREG_FN(floor,   x.w[0] += (uint32_t)fpr_floor(x.f); INREG_BAR(x); )

static void
inreg_row(const char *name, uint32_t (*fn)(void))
{
	uint32_t best = min_of(fn);
	uint32_t v = (best > inreg_cal) ? (best - inreg_cal) : 0;
	char line[64];

	sprintf(line, "%-20s %6lu  spread %lu  (16 in a row: %lu)", name,
		(unsigned long)((v + INREG_N / 2) / INREG_N),
		(unsigned long)((g_spread + INREG_N / 2) / INREG_N),
		(unsigned long)v);
	hal_send_str(line);
}

static void
bench_inreg(void)
{
	hal_send_str("--- inlined fpr operations, 16 in a row, per operation ---");
	inreg_cal = min_of(&inreg_nop);
#ifdef FNDSA_FPR_ADD_INLINE
	inreg_row("fpr_add",     &inreg_add);
#endif
#ifdef FNDSA_FPR_SUB_INLINE
	inreg_row("fpr_sub",     &inreg_sub);
#endif
#ifdef FNDSA_FPR_ADD_SUB_INLINE
	inreg_row("fpr_add_sub", &inreg_add_sub);
#endif
	inreg_row("fpr_neg",     &inreg_neg);
	inreg_row("fpr_half",    &inreg_half);
#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264
	inreg_row("fpr_of32",    &inreg_of32);
#endif
	inreg_row("fpr_rint",    &inreg_rint);
	inreg_row("fpr_floor",   &inreg_floor);
}

static void
bench_poly(unsigned logn)
{
	char hdr[64];
	size_t hn = ((size_t)1 << logn) >> 1;

	setup(logn);
	setup_fft_input(logn);
	sprintf(hdr, "--- polynomial routines, logn=%u, one call ---", logn);
	hal_send_str(hdr);

	/* Classic and compact FFT on the same input. */
	POLY("FFT old (int src)", bench_fpoly_FFT(logn, work), reseed_fft_old);
#if FNDSA_FFT_NEW

	POLY("FFT int32 compact",
		bench_fpoly_FFT_int32(logn, work, fft_int32_alias(logn)),
		reseed_fft_int32);
	POLY("FFT int32 selfadj",
		bench_fpoly_FFT_int32_selfadj(logn, work,
			fft_int32_selfadj_alias(logn)),
		reseed_fft_int32_selfadj);
	/* Including source widening, for the int8_t basis and uint16_t hm
	   inputs. */
	POLY("FFT int8 direct",
		bench_fpoly_FFT_int8(logn, work, small_src), reseed);
	POLY("FFT uint16 direct",
		bench_fpoly_FFT_uint16(logn, work, uint16_src), reseed);
	POLY("setup int8 -> int32",
		bench_expand_int8(logn, work, small_src), reseed);
	POLY("setup uint16 -> int32",
		bench_expand_uint16(logn, work, uint16_src), reseed);
	POLY("setup uint16 -> fpr",
		bench_expand_uint16_fpr(logn, work, uint16_src), reseed);
#endif

	POLY("fpoly_FFT",        bench_fpoly_FFT(logn, work), reseed);
	POLY("fpoly_iFFT",       bench_fpoly_iFFT(logn, work), reseed);
	POLY("fpoly_set_small",  bench_fpoly_set_small(logn, work, small_src),
		reseed);
	POLY("fpoly_add",        bench_fpoly_add(logn, work, operand), reseed);
	POLY("fpoly_sub",        bench_fpoly_sub(logn, work, operand), reseed);
	POLY("fpoly_neg",        bench_fpoly_neg(logn, work), reseed);
	POLY("fpoly_mul_fft",    bench_fpoly_mul_fft(logn, work, operand), reseed);
	POLY("fpoly_split_fft",  bench_fpoly_split_fft(logn, w00, w01, pristine),
		reseed);
	POLY("fpoly_merge_fft",  bench_fpoly_merge_fft(logn, work, bb00, bb01),
		reseed);

	/* The Gram matrix is built as in signing, so LDL gets real input. */
	POLY("fpoly_gram_fft",   bench_fpoly_gram_fft(logn), reseed);
	POLY("fpoly_LDL_fft",    bench_fpoly_LDL_fft(logn, w00, w01, w10),
		reseed_gram);
	POLY("fpoly_split_selfadj",
		bench_fpoly_split_selfadj(logn, w01, w10, w00), reseed_gram);

	/* mulconst and the wrapper-based ones use the empty-wrapper baseline. */
	{
		uint32_t best = 0xFFFFFFFF, worst = 0;
		for (int rep = 0; rep < REPS; rep++) {
			reseed(logn);
			uint32_t c = bench_fpoly_mulconst(logn, work);
			if (c < best) {
				best = c;
			}
			if (c > worst) {
				worst = c;
			}
		}
		report_j("fpoly_mulconst", best, cal_base0, worst - best);
	}
	(void)hn;
}

#ifdef FPRBENCH_CT
/* ---------------------------------------------------------------------- */
/* Constant-time check ('make pqm4-ctbench'). Each routine runs on edge
   cases and random inputs from its signing domain, interrupts masked.
   The M4 has no data cache, so a constant-time routine takes the same
   cycles on every input; a row passes if all counts agree. Each count is
   the minimum of CT_R runs, which drops prefetch jitter ("jit"). */

#ifndef CT_N
#define CT_N        10000
#endif
#ifndef CT_R
#define CT_R        8
#endif
#ifndef CT_POLY_N
#define CT_POLY_N   100
#endif

extern void fndsa_bw_add(void);
extern void fndsa_bw_sub(void);
extern void fndsa_bw_add_sub(void);
extern void fndsa_bw_mul(void);
extern void fndsa_bw_sqr(void);
extern void fndsa_bw_div(void);
extern void fndsa_bw_inv(void);
extern void fndsa_bw_recip_x(void);
extern void fndsa_bw_mul_recip(void);
extern void fndsa_bw_sqrt(void);
extern void fndsa_bw_sqrt_sigma(void);
extern void fndsa_bw_fpc_mul(void);
extern void fndsa_bw_neg(void);
extern void fndsa_bw_half(void);
extern void fndsa_bw_double(void);
extern void fndsa_bw_div2e(void);
extern void fndsa_bw_mul2e(void);
extern void fndsa_bw_of(void);
extern void fndsa_bw_of32(void);
extern void fndsa_bw_rint(void);
extern void fndsa_bw_floor(void);
extern void fndsa_bw_trunc(void);
#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264 && FNDSA_TWIDDLE_U64
extern void fndsa_bw_fpc_mul_gm(void);
extern void fndsa_bw_fpc_mul_gmj(void);
extern void fndsa_bw_mul_gm(void);
#endif
extern int32_t fndsa_gaussian0_helper(uint64_t lo, uint32_t hi);

static uint64_t ct_state = 0x9E3779B97F4A7C15ull;

static uint64_t
ct_rand(void)
{
	/* xorshift64*, plenty for test operands */
	ct_state ^= ct_state >> 12;
	ct_state ^= ct_state << 25;
	ct_state ^= ct_state >> 27;
	return ct_state * 0x2545F4914F6CDD1Dull;
}

/* Operand domains on the signing path. */
enum { CT_ANY, CT_POS, CT_RND };

/* The k-th operand of a domain: edge cases, then random values. */
static fpr
ct_operand(int dom, unsigned k)
{
	uint64_t r = ct_rand();
	int64_t m = (int64_t)(r >> 11);
	int sh = (int)(ct_rand() % 64);

	switch (dom) {
	case CT_POS:
		/* divisors and square-root inputs, [1, 2^29) */
		switch (k) {
		case 0: return fpr_of(1);
		case 1: return fpr_of(12289);
		case 2: return fpr_scaled(((int64_t)1 << 53) - 1, -24);
		case 3: return fpr_scaled(((int64_t)1 << 52) + 1, -52);
		}
		return fpr_scaled(m | ((int64_t)1 << 52), -52 + sh % 29);
	case CT_RND:
		/* rounding inputs, |x| < 2^29, half-integers included */
		switch (k) {
		case 0: return fpr_of(0);
		case 1: return fpr_scaled(1, -1);
		case 2: return fpr_scaled(-1, -1);
		case 3: return fpr_scaled(5, -1);
		case 4: return fpr_scaled(-5, -1);
		case 5: return fpr_scaled(1, -64);
		case 6: return fpr_scaled(-1, -64);
		}
		if (k & 1) {
			int64_t h = (int64_t)(r >> 36) | 1;
			return fpr_scaled((r & 1) ? -h : h, -1);
		}
		m >>= 24;
		return fpr_scaled((r & 1) ? -m : m, -(sh % 40));
	default:
		/* products and sums, |x| < 2^14 down to one ulp */
		switch (k) {
		case 0: return fpr_of(0);
		case 1: return fpr_of(1);
		case 2: return fpr_of(-1);
		case 3: return fpr_scaled(1, -64);
		case 4: return fpr_scaled(-1, -64);
		case 5: return fpr_scaled(((int64_t)1 << 53) - 1, -39);
		case 6: return fpr_scaled(-(((int64_t)1 << 53) - 1), -39);
		}
		return fpr_scaled((r & 1) ? -m : m, -39 - sh);
	}
}

static inline uint32_t
ct_cycles(void)
{
	return *(volatile uint32_t *)SYSCNT_ADDR;
}

/* One row's tally. */
typedef struct {
	uint32_t lo, hi, jit;
	unsigned n;
} ct_row;

static unsigned ct_rows, ct_fails;

static void
ct_row_init(ct_row *t)
{
	t->lo = 0xFFFFFFFF;
	t->hi = 0;
	t->jit = 0;
	t->n = 0;
}

static void
ct_row_add(ct_row *t, uint32_t best, uint32_t worst)
{
	if (best < t->lo) {
		t->lo = best;
	}
	if (best > t->hi) {
		t->hi = best;
	}
	if (worst - best > t->jit) {
		t->jit = worst - best;
	}
	t->n ++;
}

static void
ct_row_report(const char *name, const ct_row *t)
{
	char line[96];
	int pass = (t->lo == t->hi);

	sprintf(line, "%-22s n=%-5u min %7lu  max %7lu  jit %lu  %s",
		name, t->n, (unsigned long)t->lo, (unsigned long)t->hi,
		(unsigned long)t->jit, pass ? "PASS" : "FAIL");
	hal_send_str(line);
	ct_rows ++;
	if (!pass) {
		ct_fails ++;
	}
}

/* Time one call of a wrapper, CT_R times; best and worst. Not inlined:
   one call site, so code alignment cannot differ between rows. */
static __attribute__((noinline)) void
ct_time(void (*fn)(void), uint32_t *best, uint32_t *worst)
{
	*best = 0xFFFFFFFF;
	*worst = 0;
	for (int r = 0; r < CT_R; r ++) {
		uint32_t t0 = ct_cycles();
		fn();
		uint32_t c = ct_cycles() - t0;
		if (c < *best) {
			*best = c;
		}
		if (c > *worst) {
			*worst = c;
		}
	}
}

/* A scalar row: operands x (domain dx) and y (domain dy) per input. */
static void
ct_scalar(const char *name, void (*fn)(void), int dx, int dy)
{
	ct_row t;
	ct_row_init(&t);
	for (unsigned k = 0; k < CT_N; k ++) {
		uint32_t best, worst;
		fndsa_bw_x = ct_operand(dx, k);
		fndsa_bw_y = ct_operand(dy, k);
		fndsa_bw_i = (int32_t)ct_rand();
		ct_time(fn, &best, &worst);
		ct_row_add(&t, best, worst);
	}
	ct_row_report(name, &t);
}

/* The reciprocal-multiply consumes a reciprocal made by FPR_RECIP_X. */
static void
ct_mul_recip(void)
{
	ct_row t;
	ct_row_init(&t);
	for (unsigned k = 0; k < CT_N; k ++) {
		uint32_t best, worst;
		fndsa_bw_x = ct_operand(CT_POS, k);
		fndsa_bw_recip_x();
		fndsa_bw_x = ct_operand(CT_ANY, k);
		ct_time(fndsa_bw_mul_recip, &best, &worst);
		ct_row_add(&t, best, worst);
	}
	ct_row_report("FPR_MUL_RECIP", &t);
}

#ifdef FPR_MUL_RECIP_P63
/* The 64-bit l10 multiply: Q0.63 result, so a = x * u with |u| < 1. */
extern fpr_rcp fndsa_bw_rcp;
static volatile uint64_t ct_p63_out;

static __attribute__((noinline)) void
ct_mul_recip_p63(void)
{
	uint64_t d63;
	fpr a = fndsa_bw_y;
	FPR_MUL_RECIP_P63(d63, a, fndsa_bw_rcp);
	ct_p63_out = d63;
}

static void
ct_mul_recip_p63_row(void)
{
	ct_row t;
	ct_row_init(&t);
	for (unsigned k = 0; k < CT_N; k ++) {
		uint32_t best, worst;
		int64_t u = (int64_t)(ct_rand() >> 2);
		switch (k) {
		case 0: u = 0; break;
		case 1: u = ((int64_t)1 << 62) - 1; break;
		case 2: u = -(((int64_t)1 << 62) - 1); break;
		case 3: u = 1; break;
		default: if (k & 1) u = -u; break;
		}
		fpr x = ct_operand(CT_POS, k);
		fndsa_bw_x = x;
		fndsa_bw_recip_x();
		fndsa_bw_y = fpr_mul(x, fpr_scaled(u, -62));
		ct_time(ct_mul_recip_p63, &best, &worst);
		ct_row_add(&t, best, worst);
	}
	ct_row_report("FPR_MUL_RECIP_P63", &t);
}
#endif

/* Operations with no fprbench wrapper. */
#define NOINLINE_CT   __attribute__((noinline))
static volatile uint64_t ct_g_lo;
static volatile uint32_t ct_g_hi;
static volatile int32_t ct_g_out;

static NOINLINE_CT void
ct_gaussian0(void)
{
	ct_g_out = fndsa_gaussian0_helper(ct_g_lo, ct_g_hi);
}

static NOINLINE_CT void
ct_mul_inv_q(void)
{
	fpr x = fndsa_bw_x;
	fndsa_bw_sink_a = FPR_MUL_INV_Q(x);
}

static NOINLINE_CT void
ct_fpc_mul_self(void)
{
	fpr d_re = fndsa_bw_x, d_im = fndsa_bw_y;
	fpr b_re = fndsa_bw_y, b_im = fndsa_bw_x;
	FPC_MUL_SELF(d_re, d_im, b_re, b_im);
	fndsa_bw_sink_a = d_re;
	fndsa_bw_sink_b = d_im;
}

static NOINLINE_CT void
ct_fpc_sqr_self(void)
{
	fpr d_re = fndsa_bw_x, d_im = fndsa_bw_y;
	FPC_SQR_SELF(d_re, d_im);
	fndsa_bw_sink_a = d_re;
	fndsa_bw_sink_b = d_im;
}

static void
ct_gauss_row(void)
{
	ct_row t;
	ct_row_init(&t);
	for (unsigned k = 0; k < CT_N; k ++) {
		uint32_t best, worst;
		ct_g_lo = ct_rand();
		ct_g_hi = (uint32_t)ct_rand() & 0xFF;
		ct_time(ct_gaussian0, &best, &worst);
		ct_row_add(&t, best, worst);
	}
	ct_row_report("gaussian0_helper", &t);
}

static void
ct_scalars(void)
{
	char line[64];

	sprintf(line, "--- fpr operations, %u inputs/row ---", CT_N);
	hal_send_str(line);
	fndsa_bw_e = 9;   /* public: the degree (sqrt_sigma), a shift count */
	ct_scalar("fpr_add",        fndsa_bw_add,        CT_ANY, CT_ANY);
	ct_scalar("fpr_sub",        fndsa_bw_sub,        CT_ANY, CT_ANY);
	ct_scalar("FPR_ADD_SUB",    fndsa_bw_add_sub,    CT_ANY, CT_ANY);
	ct_scalar("fpr_mul",        fndsa_bw_mul,        CT_ANY, CT_ANY);
	ct_scalar("fpr_sqr",        fndsa_bw_sqr,        CT_ANY, CT_ANY);
	ct_scalar("FPC_MUL",        fndsa_bw_fpc_mul,    CT_ANY, CT_ANY);
	ct_scalar("FPC_MUL_SELF",   ct_fpc_mul_self,     CT_ANY, CT_ANY);
	ct_scalar("FPC_SQR_SELF",   ct_fpc_sqr_self,     CT_ANY, CT_ANY);
#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264 && FNDSA_TWIDDLE_U64
	ct_scalar("FPC_MUL_GM",     fndsa_bw_fpc_mul_gm,  CT_ANY, CT_ANY);
	ct_scalar("FPC_MUL_GMJ",    fndsa_bw_fpc_mul_gmj, CT_ANY, CT_ANY);
	ct_scalar("fpr_mul_gm",     fndsa_bw_mul_gm,      CT_ANY, CT_ANY);
#endif
	ct_scalar("FPR_MUL_INV_Q",  ct_mul_inv_q,        CT_ANY, CT_ANY);
	ct_scalar("FPR_RECIP_X",    fndsa_bw_recip_x,    CT_POS, CT_ANY);
	ct_mul_recip();
#ifdef FPR_MUL_RECIP_P63
	ct_mul_recip_p63_row();
#endif
	ct_scalar("FPR_SQRT_SIGMA", fndsa_bw_sqrt_sigma, CT_POS, CT_ANY);
	ct_scalar("fpr_sqrt",       fndsa_bw_sqrt,       CT_POS, CT_ANY);
	ct_scalar("fpr_inv",        fndsa_bw_inv,        CT_POS, CT_ANY);
	ct_scalar("fpr_div",        fndsa_bw_div,        CT_ANY, CT_POS);
	ct_scalar("fpr_neg",        fndsa_bw_neg,        CT_ANY, CT_ANY);
	ct_scalar("fpr_half",       fndsa_bw_half,       CT_ANY, CT_ANY);
	ct_scalar("fpr_double",     fndsa_bw_double,     CT_ANY, CT_ANY);
	ct_scalar("fpr_div2e",      fndsa_bw_div2e,      CT_ANY, CT_ANY);
	ct_scalar("fpr_mul2e",      fndsa_bw_mul2e,      CT_ANY, CT_ANY);
	ct_scalar("fpr_of",         fndsa_bw_of,         CT_ANY, CT_ANY);
	ct_scalar("fpr_of32",       fndsa_bw_of32,       CT_ANY, CT_ANY);
	ct_scalar("fpr_rint",       fndsa_bw_rint,       CT_RND, CT_ANY);
	ct_scalar("fpr_floor",      fndsa_bw_floor,      CT_RND, CT_ANY);
	ct_scalar("fpr_trunc",      fndsa_bw_trunc,      CT_RND, CT_ANY);
	ct_gauss_row();
}

/* Random signing-like polynomial data: small sources in [-40, 40] (the
   key polynomials), their FFTs, and an FFT-domain basis. */
static void
ct_poly_setup(unsigned logn)
{
	size_t n = (size_t)1 << logn;
	size_t sz = n * sizeof(fpr);
	fpr *dst[6] = { pristine, operand, bb00, bb01, bb10, bb11 };

	for (int p = 0; p < 6; p ++) {
		for (size_t i = 0; i < n; i ++) {
			small_src[i] = (int8_t)((int)(ct_rand() % 81) - 40);
		}
		fpoly_set_small(logn, dst[p], small_src);
		fpoly_FFT(logn, dst[p]);
	}
	for (size_t i = 0; i < n; i ++) {
		int32_t x = (int32_t)(ct_rand() % 251) - 125;
		int_src[i] = x;
		uint16_src[i] = (uint16_t)(ct_rand() % 12289);
		fft_input[i] = fpr_of32(x);
	}
	memcpy(work, pristine, sz);
	memcpy(w00, bb00, sz);
	memcpy(w01, bb01, sz);
	memcpy(w10, bb10, sz);
}

#define CT_POLY_ROW(idx, name, expr, seedfn)   do { \
		int x_ = (idx); \
		uint32_t b_ = 0xFFFFFFFF, w_ = 0; \
		for (int r_ = 0; r_ < CT_R; r_ ++) { \
			seedfn(logn); \
			uint32_t c_ = (expr); \
			if (c_ < b_) { \
				b_ = c_; \
			} \
			if (c_ > w_) { \
				w_ = c_; \
			} \
		} \
		ct_row_add(&rows[x_], b_, w_); \
		names[x_] = (name); \
	} while (0)

#define CT_POLY_MAX   24

static void
ct_polys(unsigned logn)
{
	ct_row rows[CT_POLY_MAX];
	const char *names[CT_POLY_MAX];
	int nr = 0;
	char line[64];

	sprintf(line, "--- fpoly, logn=%u, %u inputs/row ---", logn, CT_POLY_N);
	hal_send_str(line);
	for (int i = 0; i < CT_POLY_MAX; i ++) {
		ct_row_init(&rows[i]);
		names[i] = NULL;
	}
	for (unsigned k = 0; k < CT_POLY_N; k ++) {
		int i = 0;
		ct_poly_setup(logn);
		CT_POLY_ROW(i ++, "fpoly_FFT",
			bench_fpoly_FFT(logn, work), reseed);
		CT_POLY_ROW(i ++, "fpoly_iFFT",
			bench_fpoly_iFFT(logn, work), reseed);
		CT_POLY_ROW(i ++, "fpoly_set_small",
			bench_fpoly_set_small(logn, work, small_src), reseed);
		CT_POLY_ROW(i ++, "fpoly_add",
			bench_fpoly_add(logn, work, operand), reseed);
		CT_POLY_ROW(i ++, "fpoly_sub",
			bench_fpoly_sub(logn, work, operand), reseed);
		CT_POLY_ROW(i ++, "fpoly_neg",
			bench_fpoly_neg(logn, work), reseed);
		CT_POLY_ROW(i ++, "fpoly_mul_fft",
			bench_fpoly_mul_fft(logn, work, operand), reseed);
		CT_POLY_ROW(i ++, "fpoly_split_fft",
			bench_fpoly_split_fft(logn, w00, w01, pristine), reseed);
		CT_POLY_ROW(i ++, "fpoly_merge_fft",
			bench_fpoly_merge_fft(logn, work, bb00, bb01), reseed);
		CT_POLY_ROW(i ++, "fpoly_gram_fft",
			bench_fpoly_gram_fft(logn), reseed);
		CT_POLY_ROW(i ++, "fpoly_LDL_fft",
			bench_fpoly_LDL_fft(logn, w00, w01, w10), reseed_gram);
		CT_POLY_ROW(i ++, "fpoly_split_selfadj",
			bench_fpoly_split_selfadj(logn, w01, w10, w00),
			reseed_gram);
		fndsa_bw_x = ct_operand(CT_ANY, 7 + k);
		CT_POLY_ROW(i ++, "fpoly_mulconst",
			bench_fpoly_mulconst(logn, work), reseed);
#if FNDSA_FFT_NEW
		CT_POLY_ROW(i ++, "FFT int32 compact",
			bench_fpoly_FFT_int32(logn, work, fft_int32_alias(logn)),
			reseed_fft_int32);
		CT_POLY_ROW(i ++, "FFT int32 selfadj",
			bench_fpoly_FFT_int32_selfadj(logn, work,
				fft_int32_selfadj_alias(logn)),
			reseed_fft_int32_selfadj);
		CT_POLY_ROW(i ++, "FFT int8 direct",
			bench_fpoly_FFT_int8(logn, work, small_src), reseed);
		CT_POLY_ROW(i ++, "FFT uint16 direct",
			bench_fpoly_FFT_uint16(logn, work, uint16_src), reseed);
#endif
		nr = i;
	}
	for (int i = 0; i < nr; i ++) {
		ct_row_report(names[i], &rows[i]);
	}
}

static void
ct_main(void)
{
	char line[64];

	hal_send_str("constant-time check: per-row cycle range over inputs");
	/* interrupts off: a SysTick tick would land in random rows */
	__asm__ volatile ("cpsid i" ::: "memory");
	ct_scalars();
	ct_polys(9);
	ct_polys(10);
	__asm__ volatile ("cpsie i" ::: "memory");
	sprintf(line, "ct: %u rows, %u FAIL", ct_rows, ct_fails);
	hal_send_str(line);
}
#endif

int
main(void)
{
	hal_setup(CLOCK_BENCHMARK);
	enable_cyccnt();

	hal_send_str("==========================");
	hal_send_str("fpr backend micro-benchmark");
	hal_send_str(CRYPTO_ALGNAME);

#ifdef FPRBENCH_CT
	fndsa_bw_p00 = w00;
	fndsa_bw_p01 = w01;
	fndsa_bw_p10 = w10;
	fndsa_bw_p11 = bb11;
	ct_main();
	hal_send_str("#");
	hal_send_str("==========================");
	while (1) {
		/* done */
	}
#endif

	/* Calibration. */
	cal_direct = min_of(&bench_none) - 2;   /* 'bx lr' costs 2 cycles */
	cal_base0  = min_of(&bench_base0);
	cal_base1  = min_of(&bench_base1);
	cal_base1nc = min_of(&bench_base1nc);
	cal_base2  = min_of(&bench_base2);
	cal_base2nc = min_of(&bench_base2nc);
	cal_base2s2 = min_of(&bench_base2s2);
	cal_base2s2nc = min_of(&bench_base2s2nc);
	cal_addsub = min_of(&bench_base_addsub);
	cal_recip_x = min_of(&bench_base_recip_x);
	cal_mul_recip = min_of(&bench_base_mul_recip);
	cal_sqrt_sigma = min_of(&bench_base_sqrt_sigma);

	char line[96];
	sprintf(line, "cycles: minimum of %d runs, calibration subtracted;",
		REPS);
	hal_send_str(line);
	hal_send_str("spread: maximum minus minimum");
	sprintf(line, "calibration: harness=%lu wrap=%lu ld1=%lu ld2=%lu"
		" ld2nc=%lu ld2s2=%lu adds=%lu",
		(unsigned long)cal_direct, (unsigned long)cal_base0,
		(unsigned long)cal_base1, (unsigned long)cal_base2,
		(unsigned long)cal_base2nc, (unsigned long)cal_base2s2,
		(unsigned long)cal_addsub);
	hal_send_str(line);
	sprintf(line, "sizeof(fpr)=%u", (unsigned)sizeof(fpr));
	hal_send_str(line);

	/* Wrapper pointers for fpoly_gram_fft (five arguments). */
	fndsa_bw_p00 = w00;
	fndsa_bw_p01 = w01;
	fndsa_bw_p10 = w10;
	fndsa_bw_p11 = bb11;

	bench_scalar();
	bench_direct();
	bench_inreg();
	bench_poly(9);
	bench_poly(10);

	hal_send_str("#");
	hal_send_str("==========================");
	while (1) {
		/* done */
	}
	return 0;
}
