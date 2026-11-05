/*
 * Copyright 2025 NXP
 * SPDX-License-Identifier: MIT
 */

#include "triple_float.h"
#include <math.h>

const float u = 0x1p-24;
const float u_times2_plus1 = 2.0f * u + 1.0f;
const float u_times2_plus1_minus = -(2.0f * u + 1.0f);
const float u_times2_min1 = 1.0f - 2.0f * u;
const float u_times4_plus1 = 4.0f * u + 1.0f;

static inline void two_sum(const float a, const float b, float *r, float *e) {
    float aa, bb, delta_a, delta_b, s;

    s = a + b;
    aa = s - b;
    bb = s - aa;
    delta_a = a - aa;
    delta_b = b - bb;
    *e = delta_a + delta_b;
    *r = s;
}

/**
 * r+e = -a+b
 */
static inline void two_sum_neg1(const float a, const float b, float *r,
                                float *e) {
    float aa, bb, delta_a, delta_b, s;

    s = b - a;
    aa = s - b;
    bb = s - aa;
    delta_a = -a - aa;
    delta_b = b - bb;
    *e = delta_a + delta_b;
    *r = s;
}
/**
 * r+e = -a-b
 */
static inline void two_sum_neg2(const float a, const float b, float *r,
                                float *e) {
    float aa, bb, delta_a, delta_b, s;

    s = -a - b;
    aa = s + b;
    bb = s - aa;
    delta_a = -a - aa;
    delta_b = -b - bb;
    *e = delta_a + delta_b;
    *r = s;
}

static inline void fast_two_sum(const float a, const float b, float *r,
                                float *e) {
    float z;

    *r = a + b;
    z = *r - a;
    *e = b - z;
}

/**
 * r + e = -a + b
 */
static inline void fast_two_sum_neg1(const float a, const float b, float *r,
                                     float *e) {
    float z;

    *r = b - a;
    z = *r + a;
    *e = b - z;
}

/**
 * r + e = -a + -b
 */
static inline void fast_two_sum_neg2(const float a, const float b, float *r,
                                     float *e) {
    float z;

    *r = -b - a;
    z = *r + a;
    *e = -b - z;
}
/**
 * RN(c + a*b)
 */
static inline void fma_wrapper(const float a, const float b, float c,
                               float *r) {
#if !FMA_ARMV7
    *r = fmaf(a, b, c);
#else
    __asm__("vfma.f32 %0, %1, %2" : "+t"(c) : "t"(a), "t"(b));
    *r = c;
#endif
}

/**
 * RN(-c + a*b)
 */
static inline void fnms_wrapper(const float a, const float b, float c,
                                float *r) {

#if !FMA_ARMV7
    *r = fmaf(a, b, -c);
#else
    __asm__("vfnms.f32 %0, %1, %2" : "+t"(c) : "t"(a), "t"(b));
    *r = c;
#endif
}

/**
 * RN(-c + -(a*b))
 */
static inline void fnma_wrapper(const float a, const float b, float c,
                                float *r) {

#if !FMA_ARMV7
    *r = fmaf(-a, b, -c);
#else
    __asm__("vfnma.f32 %0, %1, %2" : "+t"(c) : "t"(a), "t"(b));
    *r = c;
#endif
}

/**
 * RN(c + -(a*b))
 */
static inline void fms_wrapper(const float a, const float b, float c,
                               float *r) {

#if !FMA_ARMV7
    *r = fmaf(-a, b, c);
#else
    __asm__("vfms.f32 %0, %1, %2" : "+t"(c) : "t"(a), "t"(b));
    *r = c;
#endif
}

static inline void two_prod(const float a, const float b, float *r, float *e) {
    float x = a * b;
    fnms_wrapper(a, b, x, e);
    *r = x;
}

static inline void vec_sum(const int n, const float *x, float *e) {
    float s[n];

    s[n - 1] = x[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        two_sum(x[i], s[i + 1], &s[i], &e[i + 1]);
    }
    e[0] = s[0];
}
static inline void vec_sum3(const float x0, const float x1, const float x2,
                            float *e0, float *e1, float *e2) {
    float s1;
    two_sum(x1, x2, &s1, e2);
    two_sum(x0, s1, e0, e1);
}
static inline void vec_sum3_prod(const float x0, const float x1, const float x2,
                                 float *e0, float *e1, float *e2) {
    float s1;
    fast_two_sum(x1, x2, &s1, e2);
    fast_two_sum(x0, s1, e0, e1);
}
static inline void vec_sum3_prod_min(const float x0, const float x1,
                                     const float x2, float *e0, float *e1,
                                     float *e2) {
    float s1;
    fast_two_sum_neg2(x1, x2, &s1, e2);
    fast_two_sum_neg1(x0, s1, e0, e1);
}
static inline void vec_sum4_ct(const float x0, const float x1, const float x2,
                               const float x3, float *e0, float *e1, float *e2,
                               float *e3) {
    float s1, s2;

    fast_two_sum(x2, x3, &s2, e3);
    two_sum(x1, s2, &s1, e2);
    two_sum(x0, s1, e0, e1);
}
static inline void vec_sum4_prod(const float x0, const float x1, const float x2,
                                 const float x3, float *e0, float *e1,
                                 float *e2, float *e3) {
    float s1, s2;

    fast_two_sum(x2, x3, &s2, e3);
    fast_two_sum(x1, s2, &s1, e2);
    fast_two_sum(x0, s1, e0, e1);
}
static inline void vec_sum4_prod_min(const float x0, const float x1,
                                     const float x2, const float x3, float *e0,
                                     float *e1, float *e2, float *e3) {
    float s1, s2;

    fast_two_sum_neg2(x2, x3, &s2, e3);
    fast_two_sum_neg1(x1, s2, &s1, e2);
    fast_two_sum_neg1(x0, s1, e0, e1);
}

static inline void vec_sum5_prod(const float x0, const float x1, const float x2,
                                 const float x3, const float x4, float *e0,
                                 float *e1, float *e2, float *e3, float *e4) {
    float s1, s2, s3;

    two_sum(x3, x4, &s3, e4);
    fast_two_sum(x2, s3, &s2, e3);
    fast_two_sum(x1, s2, &s1, e2);
    fast_two_sum(x0, s1, e0, e1);
}

static inline void vec_sum6(const float x0, const float x1, const float x2,
                            const float x3, const float x4, const float x5,
                            float *e0, float *e1, float *e2, float *e3,
                            float *e4, float *e5) {
    float s1, s2, s3, s4;

    two_sum(x4, x5, &s4, e5);
    two_sum(x3, s4, &s3, e4);
    two_sum(x2, s3, &s2, e3);
    two_sum(x1, s2, &s1, e2);
    two_sum(x0, s1, e0, e1);
}
static inline void vec_sum6_ct(const float x0, const float x1, const float x2,
                               const float x3, const float x4, const float x5,
                               float *e0, float *e1, float *e2, float *e3,
                               float *e4, float *e5) {
    float s1, s2, s3, s4;

    fast_two_sum(x4, x5, &s4, e5);
    two_sum(x3, s4, &s3, e4);
    two_sum(x2, s3, &s2, e3);
    two_sum(x1, s2, &s1, e2);
    two_sum(x0, s1, e0, e1);
}

#define VSEB_CT_INIT()                                                         \
    do {                                                                       \
        fast_two_sum(e0, e1, &r0, &ep_u.x);                                    \
        mask_bit = (ep_u.i | -ep_u.i) >> 31;                                   \
        mask = -mask_bit;                                                      \
        r_u.x = r0;                                                            \
        yu[idx].i = r_u.i & mask;                                              \
        idx += mask_bit;                                                       \
        ep_u.i = (ep_u.i & mask) | (r_u.i & ~mask);                            \
    } while (0);

#define VSEB_CT_STEP(k, j)                                                     \
    do {                                                                       \
        fast_two_sum(ep_u.x, e##j, &r##k, &ep_u.x);                            \
        mask_bit = (ep_u.i | -ep_u.i) >> 31;                                   \
        mask = -mask_bit;                                                      \
        r_u.x = r##k;                                                          \
        yu[idx].i = r_u.i & mask;                                              \
        idx += mask_bit;                                                       \
        ep_u.i = (ep_u.i & mask) | (r_u.i & ~mask);                            \
    } while (0);

static inline void vseb(const int n, const int m, const float *e, float *y) {
    int j;
    float ep[n];
    float ept[n];
    float r[n];

    j = 0;
    ep[0] = e[0];

    for (int i = 0; i < n - 2; i++) {
        two_sum(ep[i], e[i + 1], &r[i], &ept[i + 1]);
        if (ept[i + 1] != 0.) {
            y[j] = r[i];
            ep[i + 1] = ept[i + 1];
            j++;
        } else {
            ep[i + 1] = r[i];
        }
    }

    two_sum(ep[n - 2], e[n - 1], &y[j], &y[j + 1]);
    for (int i = j + 2; i < m; i++) {
        y[i] = 0;
    }
}

static inline void vseb_prod(const float e0, const float e1, const float e2,
                             const float e3, float *y0, float *y1) {
    float ep[2];
    float r[2];

    fast_two_sum(e0, e1, &r[0], &ep[0]);
    if (ep[0] != 0.) {
        *y0 = r[0];
        fast_two_sum(ep[0], e2, &r[1], &ep[1]);
        if (ep[1] != 0.) {
            *y1 = r[1];
        } else {
            *y1 = r[1] + e3;
        }
    } else {
        ep[0] = r[0];
        fast_two_sum(ep[0], e2, &r[1], &ep[1]);
        if (ep[1] != 0.) {
            *y0 = r[1];
            *y1 = ep[1] + e3;

        } else {
            fast_two_sum(r[1], e3, y0, y1);
        }
    }
}

static inline void vseb_prod_fast(const float e0, const float e1,
                                  const float e2, float y[2]) {
    float ep;
    float r;

    fast_two_sum(e0, e1, &r, &ep);
    if (ep != 0.) {
        y[0] = r;
        y[1] = ep + e2;
    } else {
        fast_two_sum(r, e2, &y[0], &y[1]);
    }
}
static inline void vseb_prod_fast_ct(const float e0, const float e1,
                                     const float e2, float y[2]) {
    float r;
    funion_t ep_u, r_u;
    uint32_t mask_bit;
    uint32_t mask;
    funion_t yu[3];

    uint32_t idx = 0;

    fast_two_sum(e0, e1, &r, &ep_u.x);

    mask_bit = (ep_u.i | -ep_u.i) >> 31;
    mask = -mask_bit;
    r_u.x = r;
    yu[idx].i = r_u.i & mask;
    idx += mask_bit;
    ep_u.i = (ep_u.i & mask) | (r_u.i & ~mask);

    fast_two_sum(ep_u.x, e2, &yu[idx].x, &yu[idx + 1].x);

    y[0] = yu[0].x;
    y[1] = yu[1].x;
}
static inline void vseb_prod_fast_ct_2(const float e0, const float e1,
                                       const float e2, float y[2]) {
    funion_t ep_u, r_u, ep_u2, r_u1;
    uint32_t mask;

    fast_two_sum(e0, e1, &r_u.x, &ep_u.x);
    fast_two_sum(r_u.x, e2, &r_u1.x, &ep_u2.x);

    mask = -(ep_u.i | -ep_u.i) >> 31;
    funion_t epe2;
    epe2.x = ep_u.x + e2;

    funion_t y0, y1;
    y0.i = (r_u.i & mask) | (r_u1.i & ~mask);
    y1.i = (epe2.i & ~mask) | (ep_u2.i & mask);

    y[0] = y0.x;
    y[1] = y1.x;
}

static inline void vseb_sum(const float e0, const float e1, const float e2,
                            const float e3, const float e4, const float e5,
                            float y[3]) {
    int j;
    float ep0, ep1, ep2, ep3;
    float r0, r1, r2, r3;
    float yy[6];

    j = 0;

    fast_two_sum(e0, e1, &r0, &ep0);
    if (ep0 != 0.) {
        yy[0] = r0;
        j++;
    } else {
        ep0 = r0;
    }

    fast_two_sum(ep0, e2, &r1, &ep1);
    if (ep1 != 0.) {
        yy[j] = r1;
        j++;
    } else {
        ep1 = r1;
    }

    fast_two_sum(ep1, e3, &r2, &ep2);
    if (ep2 != 0.) {
        yy[j] = r2;
        j++;
    } else {
        ep2 = r2;
    }

    fast_two_sum(ep2, e4, &r3, &ep3);
    if (ep3 != 0.) {
        yy[j] = r3;
        j++;
    } else {
        ep3 = r3;
    }

    fast_two_sum(ep3, e5, &yy[j], &yy[j + 1]);

    if (j == 0) {
        yy[2] = 0;
    }

    y[0] = yy[0];
    y[1] = yy[1];
    y[2] = yy[2];
}

static inline void vseb_sum_ct(const float e0, const float e1, const float e2,
                               const float e3, const float e4, const float e5,
                               float y[3]) {
    int idx;
    float r0, r1, r2, r3;
    funion_t yu[6];
    funion_t ep_u;
    funion_t r_u;
    uint32_t mask_bit;
    uint32_t mask;

    // there are cases where yu[2] will not be set, it then has to be 0,
    // to keep it constant time it is initialized as 0 and overwritten if
    // required
    yu[2].i = 0;
    idx = 0;

    fast_two_sum(e0, e1, &r0, &ep_u.x);
    mask_bit = (ep_u.i | -ep_u.i) >> 31;
    mask = -mask_bit;
    r_u.x = r0;
    yu[idx].i = r_u.i & mask;
    idx += mask_bit;
    ep_u.i = (ep_u.i & mask) | (r_u.i & ~mask);

    fast_two_sum(ep_u.x, e2, &r1, &ep_u.x);
    mask_bit = (ep_u.i | -ep_u.i) >> 31;
    mask = -mask_bit;
    r_u.x = r1;
    yu[idx].i = r_u.i & mask;
    idx += mask_bit;
    ep_u.i = (ep_u.i & mask) | (r_u.i & ~mask);

    fast_two_sum(ep_u.x, e3, &r2, &ep_u.x);
    mask_bit = (ep_u.i | -ep_u.i) >> 31;
    mask = -mask_bit;
    r_u.x = r2;
    yu[idx].i = r_u.i & mask;
    idx += mask_bit;
    ep_u.i = (ep_u.i & mask) | (r_u.i & ~mask);

    fast_two_sum(ep_u.x, e4, &r3, &ep_u.x);
    mask_bit = (ep_u.i | -ep_u.i) >> 31;
    mask = -mask_bit;
    r_u.x = r3;
    yu[idx].i = r_u.i & mask;
    idx += mask_bit;
    ep_u.i = (ep_u.i & mask) | (r_u.i & ~mask);

    fast_two_sum(ep_u.x, e5, &yu[idx].x, &yu[idx + 1].x);

    y[0] = yu[0].x;
    y[1] = yu[1].x;
    y[2] = yu[2].x;
}
static inline void vseb_sum_f_ct(const float e0, const float e1, const float e2,
                                 const float e3, float y[3]) {
    int idx;
    float r0, r1;
    funion_t yu[4];
    funion_t ep_u;
    funion_t r_u;
    uint32_t mask_bit;
    uint32_t mask;

    // there are cases where yu[2] will not be set, it then has to be 0,
    // to keep it constant time it is initialized as 0 and overwritten if
    // required
    yu[2].i = 0;
    idx = 0;

    fast_two_sum(e0, e1, &r0, &ep_u.x);
    mask_bit = (ep_u.i | -ep_u.i) >> 31;
    mask = -mask_bit;
    r_u.x = r0;
    yu[idx].i = r_u.i & mask;
    idx += mask_bit;
    ep_u.i = (ep_u.i & mask) | (r_u.i & ~mask);

    fast_two_sum(ep_u.x, e2, &r1, &ep_u.x);
    mask_bit = (ep_u.i | -ep_u.i) >> 31;
    mask = -mask_bit;
    r_u.x = r1;
    yu[idx].i = r_u.i & mask;
    idx += mask_bit;
    ep_u.i = (ep_u.i & mask) | (r_u.i & ~mask);

    fast_two_sum(ep_u.x, e3, &yu[idx].x, &yu[idx + 1].x);

    y[0] = yu[0].x;
    y[1] = yu[1].x;
    y[2] = yu[2].x;
}

static inline void vseb_to_tw_u64(const float e0, const float e1,
                                  const float e2, float y[3]) {
    float ep;
    float r;

    fast_two_sum(e0, e1, &r, &ep);
    if (ep != 0.) {
        y[0] = r;
        fast_two_sum(ep, e2, &y[1], &y[2]);

    } else {
        ep = r;
        fast_two_sum(r, e2, &y[0], &y[1]);
        y[2] = 0;
    }
}
tw_fpr to_tw(const float a, const float b, const float c) {
    float d0, d1;

    two_sum(a, b, &d0, &d1);

    float xv[3];
    float e[3];
    xv[0] = d0;
    xv[1] = d1;
    xv[2] = c;

    vec_sum(3, xv, e);

    tw_fpr y;
    vseb(3, 3, e, y.x);

    return y;
}

/**
 * we know a > b > c, thus we can use quicker versions of algorithms
 */
static inline tw_fpr to_tw_u64(const float a, const float b, const float c) {
    float d0, d1, e0, e1, e2;
    tw_fpr r;

    fast_two_sum(a, b, &d0, &d1);
    vec_sum3(d0, d1, c, &e0, &e1, &e2);
    vseb_to_tw_u64(e0, e1, e2, r.x);
    return r;
}

static inline void merge_ct(const tw_fpr a, const tw_fpr b, float z[6]) {

    // swap(0, 3, z);
    funion_t zi;
    zi.x = a.x[0];
    funion_t zj;
    zj.x = b.x[0];

    int32_t mask_before_shift = (zi.i & 0x7FFFFFFF) - (zj.i & 0x7FFFFFFF);

    int diff = zj.i ^ zi.i;

    // swap if z[i] > z[j]
    zj.i = zj.i ^ (diff & (mask_before_shift >> 31));
    zi.i = zi.i ^ (diff & (mask_before_shift >> 31));

    z[0] = zi.x;
    z[3] = zj.x;

    // swap(2, 5, z);
    zi.x = a.x[2];
    zj.x = b.x[2];

    mask_before_shift = (zi.i & 0x7FFFFFFF) - (zj.i & 0x7FFFFFFF);

    diff = zj.i ^ zi.i;

    // swap if z[i] > z[j]
    zj.i = zj.i ^ (diff & (mask_before_shift >> 31));
    zi.i = zi.i ^ (diff & (mask_before_shift >> 31));

    z[2] = zi.x;
    z[5] = zj.x;

    // swap(1, 4, z);
    zi.x = a.x[1];
    zj.x = b.x[1];

    mask_before_shift = (zi.i & 0x7FFFFFFF) - (zj.i & 0x7FFFFFFF);

    diff = zj.i ^ zi.i;

    // swap if z[i] > z[j]
    zj.i = zj.i ^ (diff & (mask_before_shift >> 31));
    zi.i = zi.i ^ (diff & (mask_before_shift >> 31));

    z[1] = zi.x;
    z[4] = zj.x;

    // swap(2, 3, z);
    zi.x = z[2];
    zj.x = z[3];

    mask_before_shift = (zi.i & 0x7FFFFFFF) - (zj.i & 0x7FFFFFFF);

    diff = zj.i ^ zi.i;

    // swap if z[i] > z[j]
    zj.i = zj.i ^ (diff & (mask_before_shift >> 31));
    zi.i = zi.i ^ (diff & (mask_before_shift >> 31));

    z[2] = zi.x;
    z[3] = zj.x;

    // swap(3, 4, z);
    zi.x = z[3];
    zj.x = z[4];

    mask_before_shift = (zi.i & 0x7FFFFFFF) - (zj.i & 0x7FFFFFFF);

    diff = zj.i ^ zi.i;

    // swap if z[i] > z[j]
    zj.i = zj.i ^ (diff & (mask_before_shift >> 31));
    zi.i = zi.i ^ (diff & (mask_before_shift >> 31));

    z[3] = zi.x;
    z[4] = zj.x;

    // swap(1, 2, z);
    zi.x = z[1];
    zj.x = z[2];

    mask_before_shift = (zi.i & 0x7FFFFFFF) - (zj.i & 0x7FFFFFFF);

    diff = zj.i ^ zi.i;

    // swap if z[i] > z[j]
    zj.i = zj.i ^ (diff & (mask_before_shift >> 31));
    zi.i = zi.i ^ (diff & (mask_before_shift >> 31));

    z[1] = zi.x;
    z[2] = zj.x;
}
static inline void merge_f_ct(const tw_fpr a, const float b, float z[4]) {

    funion_t zi;
    funion_t zj;

    // swap_f(0, 3, z);
    zi.x = a.x[0];
    zj.x = b;
    int32_t mask_before_shift = (zi.i & 0x7FFFFFFF) - (zj.i & 0x7FFFFFFF);

    int diff = zj.i ^ zi.i;

    // swap if z[i] > z[j]
    zj.i = zj.i ^ (diff & (mask_before_shift >> 31));
    zi.i = zi.i ^ (diff & (mask_before_shift >> 31));

    z[0] = zi.x;
    z[3] = zj.x;

    // swap_f(1, 3, z);
    zi.x = a.x[1];
    // zj.x = z[3];
    mask_before_shift = (zi.i & 0x7FFFFFFF) - (zj.i & 0x7FFFFFFF);

    diff = zj.i ^ zi.i;

    // swap if z[i] > z[j]
    zj.i = zj.i ^ (diff & (mask_before_shift >> 31));
    zi.i = zi.i ^ (diff & (mask_before_shift >> 31));

    z[1] = zi.x;
    // z[3] = zj.x;

    // swap_f(2, 3, z);
    zi.x = a.x[2];
    mask_before_shift = (zi.i & 0x7FFFFFFF) - (zj.i & 0x7FFFFFFF);

    diff = zj.i ^ zi.i;

    // swap if z[i] > z[j]
    zj.i = zj.i ^ (diff & (mask_before_shift >> 31));
    zi.i = zi.i ^ (diff & (mask_before_shift >> 31));

    z[2] = zi.x;
    z[3] = zj.x;
}

#if !FMA_ARMV7

static inline void swap_add_sub(int i, int j, float z[6], float z2[6]) {
    funion_t zi;
    zi.x = z[i];
    funion_t zj;
    zj.x = z[j];

    funion_t zi2;
    zi2.x = z2[i];
    funion_t zj2;
    zj2.x = z2[j];

    int mask = (int32_t)((zi.i & 0x7FFFFFFF) - (zj.i & 0x7FFFFFFF)) >> 31 ;

    int diff = zj.i ^ zi.i;
    int diff2 = zj2.i ^ zi2.i;

    // swap if z[i] > z[j]
    zj.i = zj.i ^ (diff & mask);
    zi.i = zi.i ^ (diff & mask);
    zj2.i = zj2.i ^ (diff2 & mask);
    zi2.i = zi2.i ^ (diff2 & mask);

    z[i] = zi.x;
    z[j] = zj.x;
    z2[i] = zi2.x;
    z2[j] = zj2.x;
}
static inline void merge_ct_add_sub(const tw_fpr a, const tw_fpr b, float z[6],
                                    float z2[6]) {

    z[0] = a.x[0];
    z[1] = a.x[1];
    z[2] = a.x[2];
    z[3] = b.x[0];
    z[4] = b.x[1];
    z[5] = b.x[2];

    z2[0] = a.x[0];
    z2[1] = a.x[1];
    z2[2] = a.x[2];
    z2[3] = -b.x[0];
    z2[4] = -b.x[1];
    z2[5] = -b.x[2];

    swap_add_sub(0, 3, z, z2);
    swap_add_sub(2, 5, z, z2);
    swap_add_sub(2, 3, z, z2);
    swap_add_sub(1, 4, z, z2);
    swap_add_sub(3, 4, z, z2);
    swap_add_sub(1, 2, z, z2);
}
#else
static inline void merge_ct_add_sub(const tw_fpr a, const tw_fpr b,
                                             float z[6], float z2[6]) {

    // swap_add_sub2(0, 3, z, z2);
    register funion_t zi;
    register funion_t zj;
    zi.x = a.x[0];
    zj.x = b.x[0];

    register funion_t zi2;
    register funion_t zj2;
    zj2.x = -b.x[0];

    register funion_t zi_n;
    register funion_t zj_n;
    register funion_t zi2_n;
    register funion_t zj2_n;


    asm volatile("bic.w %0, %4, #-2147483648\n"
        "bic.w %1, %5, #-2147483648\n"
        "subs %0, %0, %1\n"
        "sbc %1, %1\n"
        "uadd8 %0, %1, %1\n"
        "sel %0, %5, %4\n"
        "sel %1, %4, %5\n"
        "sel %2, %6, %4\n"
        "sel %3, %4, %6\n"
        : "=&r"(zi_n), "=&r"(zj_n), "=&r"(zi2_n), "=&r"(zj2_n)
        : "r"(zi.i), "r"(zj.i), "r"(zj2.i)
        : "cc");

    z[0] = zi_n.x;
    z[3] = zj_n.x;
    z2[0] = zi2_n.x;
    z2[3] = zj2_n.x;

    // swap_add_sub2(2, 5, z, z2);
    zi.x = a.x[2];
    zj.x = b.x[2];

    zj2.x = -b.x[2];

    asm volatile("bic.w %0, %4, #-2147483648\n"
        "bic.w %1, %5, #-2147483648\n"
        "subs %0, %0, %1\n"
        "sbc %1, %1\n"
        "uadd8 %0, %1, %1\n"
        "sel %0, %5, %4\n"
        "sel %1, %4, %5\n"
        "sel %2, %6, %4\n"
        "sel %3, %4, %6\n"
        : "=&r"(zi_n), "=&r"(zj_n), "=&r"(zi2_n), "=&r"(zj2_n)
        : "r"(zi.i), "r"(zj.i), "r"(zj2.i)
        : "cc");

    z[2] = zi_n.x;
    z[5] = zj_n.x;
    z2[2] = zi2_n.x;
    z2[5] = zj2_n.x;

    // swap_add_sub2(2, 3, z, z2);
    zi.x = z[2];
    zj.x = z[3];

    zi2.x = z2[2];
    zj2.x = z2[3];

    asm volatile("bic.w %0, %4, #-2147483648\n"
        "bic.w %1, %5, #-2147483648\n"
        "subs %0, %0, %1\n"
        "sbc %1, %1\n"
        "uadd8 %0, %1, %1\n"
        "sel %0, %5, %4\n"
        "sel %1, %4, %5\n"
        "sel %2, %7, %6\n"
        "sel %3, %6, %7\n"
        : "=&r"(zi_n), "=&r"(zj_n), "=&r"(zi2_n), "=&r"(zj2_n)
        : "r"(zi.i), "r"(zj.i), "r"(zi2.i), "r"(zj2.i)
        : "cc");

    z[2] = zi_n.x;
    z[3] = zj_n.x;
    z2[2] = zi2_n.x;
    z2[3] = zj2_n.x;

    // swap_add_sub2(1, 4, z, z2);
    zi.x = a.x[1];
    zj.x = b.x[1];

    zj2.x = -b.x[1];

    asm volatile("bic.w %0, %4, #-2147483648\n"
        "bic.w %1, %5, #-2147483648\n"
        "subs %0, %0, %1\n"
        "sbc %1, %1\n"
        "uadd8 %0, %1, %1\n"
        "sel %0, %5, %4\n"
        "sel %1, %4, %5\n"
        "sel %2, %6, %4\n"
        "sel %3, %4, %6\n"
        : "=&r"(zi_n), "=&r"(zj_n), "=&r"(zi2_n), "=&r"(zj2_n)
        : "r"(zi.i), "r"(zj.i), "r"(zj2.i)
        : "cc");

    z[1] = zi_n.x;
    z[4] = zj_n.x;
    z2[1] = zi2_n.x;
    z2[4] = zj2_n.x;
    // swap_add_sub2(3, 4, z, z2);
    zi.x = z[3];
    zj.x = z[4];

    zi2.x = z2[3];
    zj2.x = z2[4];

    asm volatile("bic.w %0, %4, #-2147483648\n"
        "bic.w %1, %5, #-2147483648\n"
        "subs %0, %0, %1\n"
        "sbc %1, %1\n"
        "uadd8 %0, %1, %1\n"
        "sel %0, %5, %4\n"
        "sel %1, %4, %5\n"
        "sel %2, %7, %6\n"
        "sel %3, %6, %7\n"
        : "=&r"(zi_n), "=&r"(zj_n), "=&r"(zi2_n), "=&r"(zj2_n)
        : "r"(zi.i), "r"(zj.i), "r"(zi2.i), "r"(zj2.i)
        : "cc");

    z[3] = zi_n.x;
    z[4] = zj_n.x;
    z2[3] = zi2_n.x;
    z2[4] = zj2_n.x;
    // swap_add_sub2(1, 2, z, z2);
    zi.x = z[1];
    zj.x = z[2];

    zi2.x = z2[1];
    zj2.x = z2[2];

    asm volatile("bic.w %0, %4, #-2147483648\n"
        "bic.w %1, %5, #-2147483648\n"
        "subs %0, %0, %1\n"
        "sbc %1, %1\n"
        "uadd8 %0, %1, %1\n"
        "sel %0, %5, %4\n"
        "sel %1, %4, %5\n"
        "sel %2, %7, %6\n"
        "sel %3, %6, %7\n"
        : "=&r"(zi_n), "=&r"(zj_n), "=&r"(zi2_n), "=&r"(zj2_n)
        : "r"(zi.i), "r"(zj.i), "r"(zi2.i), "r"(zj2.i)
        : "cc");

    z[1] = zi_n.x;
    z[2] = zj_n.x;
    z2[1] = zi2_n.x;
    z2[2] = zj2_n.x;
}
#endif

void merge_noloop(const tw_fpr a, const tw_fpr b, float z[6]) {
    if (__builtin_fabsf(a.x[0]) > __builtin_fabsf(b.x[0])) {
        z[0] = a.x[0];
        if (__builtin_fabsf(a.x[1]) > __builtin_fabsf(b.x[0])) {
            z[1] = a.x[1];
            if (__builtin_fabsf(a.x[2]) > __builtin_fabsf(b.x[0])) {
                z[2] = a.x[2];
                z[3] = b.x[0];
                z[4] = b.x[1];
                z[5] = b.x[2];
            } else {
                z[2] = b.x[0];
                if (__builtin_fabsf(a.x[2]) > __builtin_fabsf(b.x[1])) {
                    z[3] = a.x[2];
                    z[4] = b.x[1];
                    z[5] = b.x[2];
                } else {
                    z[3] = b.x[1];
                    z[4] = a.x[2];
                    z[5] = b.x[2];
                }
            }

        } else {
            z[1] = b.x[0];
            if (__builtin_fabsf(a.x[1]) > __builtin_fabsf(b.x[1])) {
                z[2] = a.x[1];
                if (__builtin_fabsf(a.x[2]) > __builtin_fabsf(b.x[1])) {
                    z[3] = a.x[2];
                    z[4] = b.x[1];
                    z[5] = b.x[2];
                } else {

                    z[3] = b.x[1];
                    z[4] = a.x[2];
                    z[5] = b.x[2];
                }
            } else {
                z[2] = b.x[1];
                if (__builtin_fabsf(a.x[1]) > __builtin_fabsf(b.x[2])) {
                    z[3] = a.x[1];
                    z[4] = a.x[2];
                    z[5] = b.x[2];
                } else {
                    z[3] = b.x[2];
                    z[4] = a.x[1];
                    z[5] = a.x[2];
                }
            }
        }
    } else {
        z[0] = b.x[0];
        if (__builtin_fabsf(a.x[0]) > __builtin_fabsf(b.x[1])) {
            z[1] = a.x[0];
            if (__builtin_fabsf(a.x[1]) > __builtin_fabsf(b.x[1])) {
                z[2] = a.x[1];
                if (__builtin_fabsf(a.x[2]) > __builtin_fabsf(b.x[1])) {
                    z[3] = a.x[2];
                    z[4] = b.x[1];
                    z[5] = b.x[2];
                } else {

                    z[3] = b.x[1];
                    z[4] = a.x[2];
                    z[5] = b.x[2];
                }
            } else {
                z[2] = b.x[1];
                if (__builtin_fabsf(a.x[1]) > __builtin_fabsf(b.x[2])) {
                    z[3] = a.x[1];
                    z[4] = a.x[2];
                    z[5] = b.x[2];
                } else {
                    z[3] = b.x[2];
                    z[4] = a.x[1];
                    z[5] = a.x[2];
                }
            }
        } else {
            z[1] = b.x[1];
            if (__builtin_fabsf(a.x[0]) > __builtin_fabsf(b.x[2])) {
                z[2] = a.x[0];
                if (__builtin_fabsf(a.x[1]) > __builtin_fabsf(b.x[2])) {
                    z[3] = a.x[1];
                    z[4] = a.x[2];
                    z[5] = b.x[2];
                } else {
                    z[3] = b.x[2];
                    z[4] = a.x[1];
                    z[5] = a.x[2];
                }
            } else {
                z[2] = b.x[2];
                z[3] = a.x[0];
                z[4] = a.x[1];
                z[5] = a.x[2];
            }
        }
    }
}

static inline void merge_noloop_new_add_sub(const tw_fpr a, const tw_fpr b,
                                            float z[6], float z_sub[6]) {
    if (__builtin_fabsf(a.x[0]) > __builtin_fabsf(b.x[0])) {
        z[0] = a.x[0];
        z_sub[0] = a.x[0];
        if (__builtin_fabsf(a.x[1]) > __builtin_fabsf(b.x[0])) {
            z[1] = a.x[1];
            z_sub[1] = a.x[1];
            if (__builtin_fabsf(a.x[2]) > __builtin_fabsf(b.x[0])) {
                z[2] = a.x[2];
                z[3] = b.x[0];
                z[4] = b.x[1];
                z[5] = b.x[2];
                z_sub[2] = a.x[2];
                z_sub[3] = -b.x[0];
                z_sub[4] = -b.x[1];
                z_sub[5] = -b.x[2];
            } else {
                z[2] = b.x[0];
                z_sub[2] = -b.x[0];
                if (__builtin_fabsf(a.x[2]) > __builtin_fabsf(b.x[1])) {
                    z[3] = a.x[2];
                    z[4] = b.x[1];
                    z[5] = b.x[2];
                    z_sub[3] = a.x[2];
                    z_sub[4] = -b.x[1];
                    z_sub[5] = -b.x[2];
                } else {
                    z[3] = b.x[1];
                    z[4] = a.x[2];
                    z[5] = b.x[2];
                    z_sub[3] = -b.x[1];
                    z_sub[4] = a.x[2];
                    z_sub[5] = -b.x[2];
                }
            }

        } else {
            z[1] = b.x[0];
            z_sub[1] = -b.x[0];
            if (__builtin_fabsf(a.x[1]) > __builtin_fabsf(b.x[1])) {
                z[2] = a.x[1];
                z_sub[2] = a.x[1];
                if (__builtin_fabsf(a.x[2]) > __builtin_fabsf(b.x[1])) {
                    z[3] = a.x[2];
                    z[4] = b.x[1];
                    z[5] = b.x[2];
                    z_sub[3] = a.x[2];
                    z_sub[4] = -b.x[1];
                    z_sub[5] = -b.x[2];
                } else {

                    z[3] = b.x[1];
                    z[4] = a.x[2];
                    z[5] = b.x[2];
                    z_sub[3] = -b.x[1];
                    z_sub[4] = a.x[2];
                    z_sub[5] = -b.x[2];
                }
            } else {
                z[2] = b.x[1];
                z_sub[2] = -b.x[1];
                if (__builtin_fabsf(a.x[1]) > __builtin_fabsf(b.x[2])) {
                    z[3] = a.x[1];
                    z[4] = a.x[2];
                    z[5] = b.x[2];
                    z_sub[3] = a.x[1];
                    z_sub[4] = a.x[2];
                    z_sub[5] = -b.x[2];
                } else {
                    z[3] = b.x[2];
                    z[4] = a.x[1];
                    z[5] = a.x[2];
                    z_sub[3] = -b.x[2];
                    z_sub[4] = a.x[1];
                    z_sub[5] = a.x[2];
                }
            }
        }
    } else {
        z[0] = b.x[0];
        z_sub[0] = -b.x[0];
        if (__builtin_fabsf(a.x[0]) > __builtin_fabsf(b.x[1])) {
            z[1] = a.x[0];
            z_sub[1] = a.x[0];
            if (__builtin_fabsf(a.x[1]) > __builtin_fabsf(b.x[1])) {
                z[2] = a.x[1];
                z_sub[2] = a.x[1];
                if (__builtin_fabsf(a.x[2]) > __builtin_fabsf(b.x[1])) {
                    z[3] = a.x[2];
                    z[4] = b.x[1];
                    z[5] = b.x[2];
                    z_sub[3] = a.x[2];
                    z_sub[4] = -b.x[1];
                    z_sub[5] = -b.x[2];
                } else {
                    z[3] = b.x[1];
                    z[4] = a.x[2];
                    z[5] = b.x[2];
                    z_sub[3] = -b.x[1];
                    z_sub[4] = a.x[2];
                    z_sub[5] = -b.x[2];
                }
            } else {
                z[2] = b.x[1];
                z_sub[2] = -b.x[1];
                if (__builtin_fabsf(a.x[1]) > __builtin_fabsf(b.x[2])) {
                    z[3] = a.x[1];
                    z[4] = a.x[2];
                    z[5] = b.x[2];
                    z_sub[3] = a.x[1];
                    z_sub[4] = a.x[2];
                    z_sub[5] = -b.x[2];
                } else {
                    z[3] = b.x[2];
                    z[4] = a.x[1];
                    z[5] = a.x[2];
                    z_sub[3] = -b.x[2];
                    z_sub[4] = a.x[1];
                    z_sub[5] = a.x[2];
                }
            }
        } else {
            z[1] = b.x[1];
            z_sub[1] = -b.x[1];
            if (__builtin_fabsf(a.x[0]) > __builtin_fabsf(b.x[2])) {
                z[2] = a.x[0];
                z_sub[2] = a.x[0];
                if (__builtin_fabsf(a.x[1]) > __builtin_fabsf(b.x[2])) {
                    z[3] = a.x[1];
                    z[4] = a.x[2];
                    z[5] = b.x[2];
                    z_sub[3] = a.x[1];
                    z_sub[4] = a.x[2];
                    z_sub[5] = -b.x[2];
                } else {
                    z[3] = b.x[2];
                    z[4] = a.x[1];
                    z[5] = a.x[2];
                    z_sub[3] = -b.x[2];
                    z_sub[4] = a.x[1];
                    z_sub[5] = a.x[2];
                }
            } else {
                z[2] = b.x[2];
                z[3] = a.x[0];
                z[4] = a.x[1];
                z[5] = a.x[2];
                z_sub[2] = -b.x[2];
                z_sub[3] = a.x[0];
                z_sub[4] = a.x[1];
                z_sub[5] = a.x[2];
            }
        }
    }
}

inline tw_fpr tw_sum(const tw_fpr a, const tw_fpr b) {
    float e0, e1, e2, e3, e4, e5, z[6];
    tw_fpr r;

    merge_noloop(a, b, z);
    vec_sum6(z[0], z[1], z[2], z[3], z[4], z[5], &e0, &e1, &e2, &e3, &e4, &e5);

    vseb_sum(e0, e1, e2, e3, e4, e5, r.x);
    return r;
}

tw_fpr tw_sum_ct(const tw_fpr a, const tw_fpr b) {
    float e0, e1, e2, e3, e4, e5, z[6];
    tw_fpr r;

    merge_ct(a, b, z);
    vec_sum6_ct(z[0], z[1], z[2], z[3], z[4], z[5], &e0, &e1, &e2, &e3, &e4,
                &e5);

    vseb_sum_ct(e0, e1, e2, e3, e4, e5, r.x);
    return r;
}

tw_fpr tw_sum_f_ct(const tw_fpr a, const float b) {
    float e0, e1, e2, e3, z[4];
    tw_fpr r;

    merge_f_ct(a, b, z);
    vec_sum4_ct(z[0], z[1], z[2], z[3], &e0, &e1, &e2, &e3);

    vseb_sum_f_ct(e0, e1, e2, e3, r.x);
    return r;
}

void tw_add_sub(const tw_fpr a, const tw_fpr b, tw_fpr *ra, tw_fpr *rb) {
    float e0, e1, e2, e3, e4, e5, e0_min, e1_min, e2_min, e3_min, e4_min,
        e5_min, z[6], z_min[6];

    merge_noloop_new_add_sub(a, b, z, z_min);
    vec_sum6(z[0], z[1], z[2], z[3], z[4], z[5], &e0, &e1, &e2, &e3, &e4, &e5);

    vec_sum6(z_min[0], z_min[1], z_min[2], z_min[3], z_min[4], z_min[5],
             &e0_min, &e1_min, &e2_min, &e3_min, &e4_min, &e5_min);

    vseb_sum(e0, e1, e2, e3, e4, e5, ra->x);
    vseb_sum(e0_min, e1_min, e2_min, e3_min, e4_min, e5_min, rb->x);
}

void tw_add_sub_ct(const tw_fpr a, const tw_fpr b, tw_fpr *ra, tw_fpr *rb) {
    float e0, e1, e2, e3, e4, e5, e0_min, e1_min, e2_min, e3_min, e4_min,
        e5_min, z[6], z_min[6];

    merge_ct_add_sub(a, b, z, z_min);
    vec_sum6(z[0], z[1], z[2], z[3], z[4], z[5], &e0, &e1, &e2, &e3, &e4, &e5);
    vec_sum6(z_min[0], z_min[1], z_min[2], z_min[3], z_min[4], z_min[5],
             &e0_min, &e1_min, &e2_min, &e3_min, &e4_min, &e5_min);

    vseb_sum_ct(e0, e1, e2, e3, e4, e5, ra->x);
    vseb_sum_ct(e0_min, e1_min, e2_min, e3_min, e4_min, e5_min, rb->x);
}

inline tw_fpr tw_prod_acc(const tw_fpr a, const tw_fpr b) {
    float z00_plus, z00_min, z01_plus, z01_min, z10_plus, z10_min, b0, b1, b2,
        c, z31, z32, z3, e0, e1, e2, e3, e4;
    tw_fpr r;

    two_prod(a.x[0], b.x[0], &z00_plus, &z00_min);
    two_prod(a.x[0], b.x[1], &z01_plus, &z01_min);
    two_prod(a.x[1], b.x[0], &z10_plus, &z10_min);

    vec_sum3(z00_min, z01_plus, z10_plus, &b0, &b1, &b2);

    fma_wrapper(a.x[1], b.x[1], b2, &c);
    fma_wrapper(a.x[0], b.x[2], z10_min, &z31);
    fma_wrapper(a.x[2], b.x[0], z01_min, &z32);
    z3 = z31 + z32;

    vec_sum5_prod(z00_plus, b0, b1, c, z3, &e0, &e1, &e2, &e3, &e4);

    r.x[0] = e0;

    vseb_prod(e1, e2, e3, e4, &r.x[1], &r.x[2]);
    return r;
}

inline tw_fpr tw_prod_fast(const tw_fpr a, const tw_fpr b) {
    float z00_plus, z00_min, z01_plus, z01_min, z10_plus, z10_min, b0, b1, b2,
        c, z31, z32, z3, e0, e1, e2, e3, s3;
    tw_fpr r;

    two_prod(a.x[0], b.x[0], &z00_plus, &z00_min);
    two_prod(a.x[0], b.x[1], &z01_plus, &z01_min);
    two_prod(a.x[1], b.x[0], &z10_plus, &z10_min);

    vec_sum3(z00_min, z01_plus, z10_plus, &b0, &b1, &b2);

    fma_wrapper(a.x[1], b.x[1], b2, &c);
    fma_wrapper(a.x[0], b.x[2], z10_min, &z31);
    fma_wrapper(a.x[2], b.x[0], z01_min, &z32);
    z3 = z31 + z32;
    s3 = c + z3;

    vec_sum4_prod(z00_plus, b0, b1, s3, &e0, &e1, &e2, &e3);

    r.x[0] = e0;

    vseb_prod_fast(e1, e2, e3, &r.x[1]);
    return r;
}

inline tw_fpr tw_prod_fast_ct(const tw_fpr a, const tw_fpr b) {
    float z00_plus, z00_min, z01_plus, z01_min, z10_plus, z10_min, b0, b1, b2,
        c, z31, z32, z3, e0, e1, e2, e3, s3;
    tw_fpr r;

    two_prod(a.x[0], b.x[0], &z00_plus, &z00_min);
    two_prod(a.x[0], b.x[1], &z01_plus, &z01_min);
    two_prod(a.x[1], b.x[0], &z10_plus, &z10_min);

    vec_sum3(z00_min, z01_plus, z10_plus, &b0, &b1, &b2);

    fma_wrapper(a.x[1], b.x[1], b2, &c);
    fma_wrapper(a.x[0], b.x[2], z10_min, &z31);
    fma_wrapper(a.x[2], b.x[0], z01_min, &z32);
    z3 = z31 + z32;
    s3 = c + z3;

    vec_sum4_prod(z00_plus, b0, b1, s3, &e0, &e1, &e2, &e3);

    r.x[0] = e0;

    vseb_prod_fast_ct_2(e1, e2, e3, &r.x[1]);
    return r;
}

static inline tw_fpr tw_prod_acc2(const float a0, const float a1,
                                  const tw_fpr b) {
    float z00_plus, z00_min, z01_plus, z01_min, z10_plus, z10_min, b0, b1, b2,
        c, z31, z3, e0, e1, e2, e3, e4;
    tw_fpr r;

    two_prod(a0, b.x[0], &z00_plus, &z00_min);
    two_prod(a0, b.x[1], &z01_plus, &z01_min);
    two_prod(a1, b.x[0], &z10_plus, &z10_min);

    vec_sum3(z00_min, z01_plus, z10_plus, &b0, &b1, &b2);

    fma_wrapper(a1, b.x[1], b2, &c);
    fma_wrapper(a0, b.x[2], z10_min, &z31);
    z3 = z31 + z01_min;

    vec_sum5_prod(z00_plus, b0, b1, c, z3, &e0, &e1, &e2, &e3, &e4);
    r.x[0] = e0;
    vseb_prod(e1, e2, e3, e4, &r.x[1], &r.x[2]);
    return r;
}
static inline tw_fpr tw_prod_fast_2(const float a0, const float a1,
                                    const tw_fpr b) {
    float z00_plus, z00_min, z01_plus, z01_min, z10_plus, z10_min, b0, b1, b2,
        c, z31, z3, e0, e1, e2, e3, s3;
    tw_fpr r;

    two_prod(a0, b.x[0], &z00_plus, &z00_min);
    two_prod(a0, b.x[1], &z01_plus, &z01_min);
    two_prod(a1, b.x[0], &z10_plus, &z10_min);

    vec_sum3(z00_min, z01_plus, z10_plus, &b0, &b1, &b2);

    fma_wrapper(a1, b.x[1], b2, &c);
    fma_wrapper(a0, b.x[2], z10_min, &z31);
    z3 = z31 + z01_min;
    s3 = c + z3;

    vec_sum4_prod(z00_plus, b0, b1, s3, &e0, &e1, &e2, &e3);
    r.x[0] = e0;
    vseb_prod_fast(e1, e2, e3, &r.x[1]);
    return r;
}

static inline tw_fpr tw_prod_fast_2_ct(const float a0, const float a1,
                                       const tw_fpr b) {
    float z00_plus, z00_min, z01_plus, z01_min, z10_plus, z10_min, b0, b1, b2,
        c, z31, z3, e0, e1, e2, e3, s3;
    tw_fpr r;

    two_prod(a0, b.x[0], &z00_plus, &z00_min);
    two_prod(a0, b.x[1], &z01_plus, &z01_min);
    two_prod(a1, b.x[0], &z10_plus, &z10_min);

    vec_sum3(z00_min, z01_plus, z10_plus, &b0, &b1, &b2);

    fma_wrapper(a1, b.x[1], b2, &c);
    fma_wrapper(a0, b.x[2], z10_min, &z31);
    z3 = z31 + z01_min;
    s3 = c + z3;

    vec_sum4_prod(z00_plus, b0, b1, s3, &e0, &e1, &e2, &e3);
    r.x[0] = e0;
    vseb_prod_fast_ct(e1, e2, e3, &r.x[1]);
    return r;
}
tw_fpr tw_prod_fast_f_ct(const float a0, const tw_fpr b) {
    float z00_plus, z00_min, z01_plus, z01_min, b0, b1, z3, e0, e1, e2, e3;
    tw_fpr r;

    two_prod(a0, b.x[0], &z00_plus, &z00_min);
    two_prod(a0, b.x[1], &z01_plus, &z01_min);

    two_sum(z00_min, z01_plus, &b0, &b1);

    z3 = a0 * b.x[2];
    z3 += z01_min;

    vec_sum4_prod(z00_plus, b0, b1, z3, &e0, &e1, &e2, &e3);
    r.x[0] = e0;
    vseb_prod_fast_ct(e1, e2, e3, &r.x[1]);
    return r;
}
static inline void tw_2_min_prod_acc(const float b0, const float b1,
                                     const tw_fpr x, float i[2]) {
    float z00_plus, z00_min, z01_plus, z01_min, z10_plus, z10_min, b0t, b1t,
        b2t, c, z31, z3, e1, e2, e3, e4, s1, r1;

    two_prod(b0, x.x[0], &z00_plus, &z00_min);
    two_prod(b0, x.x[1], &z01_plus, &z01_min);
    two_prod(b1, x.x[0], &z10_plus, &z10_min);

    vec_sum3(z00_min, z01_plus, z10_plus, &b0t, &b1t, &b2t);

    fma_wrapper(b1, x.x[1], b2t, &c);
    fma_wrapper(b0, x.x[2], z10_min, &z31);
    // fault in algorithm specification, says z32 instead of z31 here
    z3 = z31 + z01_min;

    vec_sum4_prod_min(b0t, b1t, c, z3, &s1, &e2, &e3, &e4);
    fast_two_sum_neg1(z00_plus, s1, &r1, &e1);

    vseb_prod(e1, e2, e3, e4, &i[0], &i[1]);
}

static inline void tw_2_min_prod_fast(const float b0, const float b1,
                                      const tw_fpr x, float i[2]) {
    float z00_plus, z00_min, z01_plus, z01_min, z10_plus, z10_min, b0t, b1t,
        b2t, c, z31, z3, e1, e2, e3, s1, r1;

    two_prod(b0, x.x[0], &z00_plus, &z00_min);
    two_prod(b0, x.x[1], &z01_plus, &z01_min);
    two_prod(b1, x.x[0], &z10_plus, &z10_min);

    vec_sum3(z00_min, z01_plus, z10_plus, &b0t, &b1t, &b2t);

    fma_wrapper(b1, x.x[1], b2t, &c);
    fma_wrapper(b0, x.x[2], z10_min, &z31);
    // fault in algorithm specification, says z32 instead of z31 here
    z3 = z31 + z01_min;
    float s3 = c + z3;

    vec_sum3_prod_min(b0t, b1t, s3, &s1, &e2, &e3);
    fast_two_sum_neg1(z00_plus, s1, &r1, &e1);

    vseb_prod_fast(e1, e2, e3, i);
}
static inline void tw_2_min_prod_fast_ct(const float b0, const float b1,
                                         const tw_fpr x, float i[2]) {
    float z00_plus, z00_min, z01_plus, z01_min, z10_plus, z10_min, b0t, b1t,
        b2t, c, z31, z3, e1, e2, e3, s1, r1;

    two_prod(b0, x.x[0], &z00_plus, &z00_min);
    two_prod(b0, x.x[1], &z01_plus, &z01_min);
    two_prod(b1, x.x[0], &z10_plus, &z10_min);

    vec_sum3(z00_min, z01_plus, z10_plus, &b0t, &b1t, &b2t);

    fma_wrapper(b1, x.x[1], b2t, &c);
    fma_wrapper(b0, x.x[2], z10_min, &z31);
    // fault in algorithm specification, says z32 instead of z31 here
    z3 = z31 + z01_min;
    float s3 = c + z3;

    vec_sum3_prod_min(b0t, b1t, s3, &s1, &e2, &e3);
    fast_two_sum_neg1(z00_plus, s1, &r1, &e1);

    vseb_prod_fast_ct_2(e1, e2, e3, i);
}

static inline void tw_1_point_5_min_prod_acc(const float b0, const float b1,
                                             const tw_fpr x, float i[2]) {
    float z00_plus, z00_min, z01_plus, z01_min, z10_plus, z10_min, b0t, b1t,
        b2t, c, z31, z3, e1, e2, e3, e4, s1, r1;

    two_prod(b0, x.x[0], &z00_plus, &z00_min);
    two_prod(b0, x.x[1], &z01_plus, &z01_min);
    two_prod(b1, x.x[0], &z10_plus, &z10_min);

    vec_sum3(z00_min, z01_plus, z10_plus, &b0t, &b1t, &b2t);

    fma_wrapper(b1, x.x[1], b2t, &c);
    fma_wrapper(b0, x.x[2], z10_min, &z31);
    // fault in algorithm specification, says z32 instead of z31 here
    z3 = z31 + z01_min;

    vec_sum4_prod_min(b0t, b1t, c, z3, &s1, &e2, &e3, &e4);
    // Algorithms specifies * 0.5 somewhere, but For some reason the only way it
    // is correct is by just doing it the normal way
    fast_two_sum_neg1(z00_plus, s1, &r1, &e1);

    vseb_prod(e1, e2, e3, e4, &i[0], &i[1]);
}

static inline void tw_1_point_5_min_prod_fast(const float b0, const float b1,
                                              const tw_fpr x, float i[2]) {
    float z00_plus, z00_min, z01_plus, z01_min, z10_plus, z10_min, b0t, b1t,
        b2t, c, z31, z3, e1, e2, e3, s1, r1, s3;

    two_prod(b0, x.x[0], &z00_plus, &z00_min);
    two_prod(b0, x.x[1], &z01_plus, &z01_min);
    two_prod(b1, x.x[0], &z10_plus, &z10_min);

    vec_sum3(z00_min, z01_plus, z10_plus, &b0t, &b1t, &b2t);

    fma_wrapper(b1, x.x[1], b2t, &c);
    fma_wrapper(b0, x.x[2], z10_min, &z31);
    // fault in algorithm specification, says z32 instead of z31 here
    z3 = z31 + z01_min;
    s3 = c + z3;

    vec_sum3_prod_min(b0t, b1t, s3, &s1, &e2, &e3);
    fast_two_sum_neg1(z00_plus, s1, &r1, &e1);

    vseb_prod_fast(e1, e2, e3, i);
}
static inline void tw_1_point_5_min_prod_fast_ct(const float b0, const float b1,
                                                 const tw_fpr x, float i[2]) {
    float z00_plus, z00_min, z01_plus, z01_min, z10_plus, z10_min, b0t, b1t,
        b2t, c, z31, z3, e1, e2, e3, s1, r1, s3;

    two_prod(b0, x.x[0], &z00_plus, &z00_min);
    two_prod(b0, x.x[1], &z01_plus, &z01_min);
    two_prod(b1, x.x[0], &z10_plus, &z10_min);

    vec_sum3(z00_min, z01_plus, z10_plus, &b0t, &b1t, &b2t);

    fma_wrapper(b1, x.x[1], b2t, &c);
    fma_wrapper(b0, x.x[2], z10_min, &z31);
    // fault in algorithm specification, says z32 instead of z31 here
    z3 = z31 + z01_min;
    s3 = c + z3;

    vec_sum3_prod_min(b0t, b1t, s3, &s1, &e2, &e3);
    fast_two_sum_neg1(z00_plus, s1, &r1, &e1);

    vseb_prod_fast_ct(e1, e2, e3, i);
}

static inline tw_fpr tw_prod_fast_reci_y(const float b0, const float b1,
                                         const float i1, const float i2) {
    tw_fpr r;
    float z01_plus, z01_min, b0t, b1t, z31, z3, s3, e0, e1, e2;

    two_prod(b0, i1, &z01_plus, &z01_min);
    fast_two_sum(b1, z01_plus, &b0t, &b1t);

    fma_wrapper(b1, i1, z01_min, &z31);
    fma_wrapper(b0, i2, z31, &z3);
    s3 = b1t + z3;

    vec_sum3_prod(b0, b0t, s3, &e0, &e1, &e2);
    r.x[0] = e0;

    fast_two_sum(e1, e2, &r.x[1], &r.x[2]);
    return r;
}

static inline tw_fpr tw_prod_fast_div_y(const tw_fpr a, const float i[2]) {
    float z01_plus, z01_min, b0t, b1t, z31, z32, z3, s3, e0, e1, e2;
    tw_fpr r;

    two_prod(a.x[0], i[0], &z01_plus, &z01_min);
    fast_two_sum(a.x[1], z01_plus, &b0t, &b1t);

    fma_wrapper(a.x[1], i[0], z01_min, &z31);
    fma_wrapper(a.x[0], i[1], z31, &z32);
    z3 = z32 + b1t;
    s3 = z3 + a.x[2];

    vec_sum3_prod(a.x[0], b0t, s3, &e0, &e1, &e2);
    r.x[0] = e0;

    fast_two_sum(e1, e2, &r.x[1], &r.x[2]);
    return r;
}

tw_fpr tw_reci(const tw_fpr x) {
    float a, h11, h1, b01, b11, b12, b_hat0, b_hat1, i_hat[2];

    a = u_times2_plus1 / x.x[0];

    fma_wrapper(a, x.x[0], -u_times2_plus1, &h11);
    fma_wrapper(-a, x.x[1], -h11, &h1);
    two_prod(a, u_times2_min1, &b01, &b11);
    fma_wrapper(a, h1, b11, &b12);
    fast_two_sum(b01, b12, &b_hat0, &b_hat1);

    tw_2_min_prod_acc(b_hat0, b_hat1, x, i_hat);
    return tw_prod_fast_reci_y(b_hat0, b_hat1, i_hat[0], i_hat[1]);
}

tw_fpr tw_reci_fast(const tw_fpr x) {
    float a, h11, h1, b01, b11, b12, b_hat0, b_hat1, i_hat[2];

    a = u_times2_plus1 / x.x[0];

    fma_wrapper(a, x.x[0], -u_times2_plus1, &h11);
    fma_wrapper(-a, x.x[1], -h11, &h1);
    two_prod(a, u_times2_min1, &b01, &b11);
    fma_wrapper(a, h1, b11, &b12);
    fast_two_sum(b01, b12, &b_hat0, &b_hat1);

    tw_2_min_prod_fast(b_hat0, b_hat1, x, i_hat);
    return tw_prod_fast_reci_y(b_hat0, b_hat1, i_hat[0], i_hat[1]);
}
tw_fpr tw_reci_fast_ct(const tw_fpr x) {
    float a, h11, h1, b01, b11, b12, b_hat0, b_hat1, i_hat[2];

    a = u_times2_plus1 / x.x[0];

    fma_wrapper(a, x.x[0], -u_times2_plus1, &h11);
    fma_wrapper(-a, x.x[1], -h11, &h1);
    two_prod(a, u_times2_min1, &b01, &b11);
    fma_wrapper(a, h1, b11, &b12);
    fast_two_sum(b01, b12, &b_hat0, &b_hat1);

    tw_2_min_prod_fast_ct(b_hat0, b_hat1, x, i_hat);
    return tw_prod_fast_reci_y(b_hat0, b_hat1, i_hat[0], i_hat[1]);
}

tw_fpr tw_div(const tw_fpr z, const tw_fpr x) {
    float a, h11, h1, b01, b11, b12, b_hat0, b_hat1, i_hat[2];
    tw_fpr a_hat;

    a = u_times2_plus1 / x.x[0];

    fnms_wrapper(a, x.x[0], u_times2_plus1, &h11);
    fms_wrapper(a, x.x[1], -h11, &h1);
    two_prod(a, u_times2_min1, &b01, &b11);
    fma_wrapper(a, h1, b11, &b12);

    fast_two_sum(b01, b12, &b_hat0, &b_hat1);

    tw_2_min_prod_acc(b_hat0, b_hat1, x, i_hat);

    a_hat = tw_prod_acc2(b_hat0, b_hat1, z);
    return tw_prod_fast_div_y(a_hat, i_hat);
}

tw_fpr tw_div_fast(const tw_fpr z, const tw_fpr x) {
    float a, h11, h1, b01, b11, b12, b_hat0, b_hat1, i_hat[2];
    tw_fpr a_hat;

    a = u_times2_plus1 / x.x[0];

    fnms_wrapper(a, x.x[0], u_times2_plus1, &h11);
    fms_wrapper(a, x.x[1], -h11, &h1);
    two_prod(a, u_times2_min1, &b01, &b11);
    fma_wrapper(a, h1, b11, &b12);

    fast_two_sum(b01, b12, &b_hat0, &b_hat1);

    tw_2_min_prod_fast(b_hat0, b_hat1, x, i_hat);

    a_hat = tw_prod_fast_2(b_hat0, b_hat1, z);
    return tw_prod_fast_div_y(a_hat, i_hat);
}
tw_fpr tw_div_fast_ct(const tw_fpr z, const tw_fpr x) {
    float a, h11, h1, b01, b11, b12, b_hat0, b_hat1, i_hat[2];
    tw_fpr a_hat;

    a = u_times2_plus1 / x.x[0];

    fnms_wrapper(a, x.x[0], u_times2_plus1, &h11);
    fms_wrapper(a, x.x[1], -h11, &h1);
    two_prod(a, u_times2_min1, &b01, &b11);
    fma_wrapper(a, h1, b11, &b12);

    fast_two_sum(b01, b12, &b_hat0, &b_hat1);

    tw_2_min_prod_fast_ct(b_hat0, b_hat1, x, i_hat);

    a_hat = tw_prod_fast_2_ct(b_hat0, b_hat1, z);
    return tw_prod_fast_div_y(a_hat, i_hat);
}

tw_fpr tw_sqrt(const tw_fpr x) {
    float x0_sqrt, a, a_prime, h0_1, h11_1, h1_1, h01_2, h11_2, h0_2, h1_2, b01,
        b11, b12, b_hat0, b_hat1, b_hat_prime0, b_hat_prime1, i_hat2[2];
    tw_fpr i1;

    x0_sqrt = __builtin_sqrtf(x.x[0]);
    a = u_times4_plus1 / x0_sqrt;
    a_prime = 0.5f * a;

    two_prod(a, x.x[0], &h0_1, &h11_1);
    fma_wrapper(a, x.x[1], h11_1, &h1_1);

    two_prod(a_prime, h0_1, &h01_2, &h11_2);
    h0_2 = 1.5f - h01_2;
    fma_wrapper(a_prime, h1_1, h11_2, &h1_2);
    h1_2 = -h1_2;

    two_prod(a, h0_2, &b01, &b11);
    fma_wrapper(a, h1_2, b11, &b12);

    fast_two_sum(b01, b12, &b_hat0, &b_hat1);
    b_hat_prime0 = 0.5f * b_hat0;
    b_hat_prime1 = 0.5f * b_hat1;

    i1 = tw_prod_acc2(b_hat0, b_hat1, x);

    tw_1_point_5_min_prod_acc(b_hat_prime0, b_hat_prime1, i1, i_hat2);

    return tw_prod_fast_div_y(i1, i_hat2);
}

tw_fpr tw_sqrt_fast(const tw_fpr x) {
    float x0_sqrt, a, a_prime, h0_1, h11_1, h1_1, h01_2, h11_2, h0_2, h1_2, b01,
        b11, b12, b_hat0, b_hat1, b_hat_prime0, b_hat_prime1, i_hat2[2];
    tw_fpr i1;

    x0_sqrt = __builtin_sqrtf(x.x[0]);
    a = u_times4_plus1 / x0_sqrt;
    a_prime = 0.5f * a;

    two_prod(a, x.x[0], &h0_1, &h11_1);
    fma_wrapper(a, x.x[1], h11_1, &h1_1);

    two_prod(a_prime, h0_1, &h01_2, &h11_2);
    h0_2 = 1.5f - h01_2;
    fma_wrapper(a_prime, h1_1, h11_2, &h1_2);
    h1_2 = -h1_2;

    two_prod(a, h0_2, &b01, &b11);
    fma_wrapper(a, h1_2, b11, &b12);

    fast_two_sum(b01, b12, &b_hat0, &b_hat1);
    b_hat_prime0 = 0.5f * b_hat0;
    b_hat_prime1 = 0.5f * b_hat1;

    i1 = tw_prod_fast_2(b_hat0, b_hat1, x);

    tw_1_point_5_min_prod_fast(b_hat_prime0, b_hat_prime1, i1, i_hat2);

    return tw_prod_fast_div_y(i1, i_hat2);
}
tw_fpr tw_sqrt_fast_ct(const tw_fpr x) {
    float x0_sqrt, a, a_prime, h0_1, h11_1, h1_1, h01_2, h11_2, h0_2, h1_2, b01,
        b11, b12, b_hat0, b_hat1, b_hat_prime0, b_hat_prime1, i_hat2[2];
    tw_fpr i1;

    x0_sqrt = __builtin_sqrtf(x.x[0]);
    a = u_times4_plus1 / x0_sqrt;
    a_prime = 0.5f * a;

    two_prod(a, x.x[0], &h0_1, &h11_1);
    fma_wrapper(a, x.x[1], h11_1, &h1_1);

    two_prod(a_prime, h0_1, &h01_2, &h11_2);
    h0_2 = 1.5f - h01_2;
    fma_wrapper(a_prime, h1_1, h11_2, &h1_2);
    h1_2 = -h1_2;

    two_prod(a, h0_2, &b01, &b11);
    fma_wrapper(a, h1_2, b11, &b12);

    fast_two_sum(b01, b12, &b_hat0, &b_hat1);
    b_hat_prime0 = 0.5f * b_hat0;
    b_hat_prime1 = 0.5f * b_hat1;

    i1 = tw_prod_fast_2_ct(b_hat0, b_hat1, x);

    tw_1_point_5_min_prod_fast_ct(b_hat_prime0, b_hat_prime1, i1, i_hat2);

    return tw_prod_fast_div_y(i1, i_hat2);
}

static inline tw_fpr form_from_u64(uint64_t x) {

    tw_fpr r;
    funion_t xf;
    uint32_t s = (x & 0x8000000000000000) >> 32;
    // set sign
    // set exponent by getting e of double precision, - 1024 to get e from
    // biased double preicion , + 127 to get biased e for single precision
    uint32_t e = (((x >> 52) & 0x7FF) - 1023 + 127);
    xf.i = s ^ e << 23;

    // set precision
    xf.i ^= (x >> 29) & 0x7FFFFF;

    // uint64_t x_tmp = x_tmp << (12 + 23);

    // get the last 29 precision bits, and shift them to the left;
    uint32_t p_leftover = (x & 0x1FFFFFFF) << 3;
    int b;

    funion_t xf2;
    if (p_leftover != 0) {
        // computes the amount of leading 0's untill the first 1, starting at
        // msb
        b = __builtin_clz(p_leftover);

        // set s
        xf2.i = s;

        // old exponent - 23, because 23 p's already used in xf.
        // -b because those were 0
        e = e - 23 - b;
        // -1 because of implicit 1?
        e = e - 1;

        // if (e < 0) {
        //     xf2.x = 0;
        // } else {
        xf2.i ^= e << 23;

        // shift untill after the first 1, as that one is implicit in
        // binary32
        p_leftover = p_leftover << (b + 1);

        // first 9 bits are sign and e, then the p values should come
        xf2.i ^= p_leftover >> 9;

        // next 23 bits used in xf2
        p_leftover <<= 23;
        // }
    } else {
        xf2.x = 0.0;
    }

    funion_t xf3;

    if (p_leftover != 0) {
        b = __builtin_clz(p_leftover);

        // old exponent - 23, because 23 p's already used in xf.
        // -b because those were 0
        e = e - 23 - b;
        // -1 because of implicit 1?
        e = e - 1;

        // if (e > 255) {
        //     xf3.x = 0;
        // } else {
        //     // set s
        xf3.i = s;
        // set e
        xf3.i ^= e << 23;

        p_leftover <<= (b + 1);

        // set p
        xf3.i ^= p_leftover >> 9;
        // }
    } else {
        xf3.x = 0.;
    }
    if (isnanf(xf.x) || isinff(xf3.x)) {
        xf.x = 0.0f;
    }
    if (isnanf(xf2.x) || isinff(xf3.x)) {
        xf2.x = 0.0f;
    }
    if (isnanf(xf3.x) || isinff(xf3.x)) {
        xf3.x = 0.0f;
    }

    r.x[0] = xf.x;
    r.x[1] = xf2.x;
    r.x[2] = xf3.x;
    return r;
}

tw_fpr create_from_u64(uint64_t x) {
    tw_fpr ri = form_from_u64(x);

    return to_tw_u64(ri.x[0], ri.x[1], ri.x[2]);
}

tw_fpr tw_sub(tw_fpr x, tw_fpr y) {
    y.x[0] = -y.x[0];
    y.x[1] = -y.x[1];
    y.x[2] = -y.x[2];
    return tw_sum(x, y);
}

tw_fpr tw_sub_ct(tw_fpr x, tw_fpr y) {
    y.x[0] = -y.x[0];
    y.x[1] = -y.x[1];
    y.x[2] = -y.x[2];
    return tw_sum_ct(x, y);
}
