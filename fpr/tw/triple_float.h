#ifndef FNDSA_FPR_TW_TRIPLE_FLOAT_H__
#define FNDSA_FPR_TW_TRIPLE_FLOAT_H__

/*
 * Copyright 2025 NXP
 * SPDX-License-Identifier: MIT
 */

#include <stdint.h>

#define FMA_ARMV7 0

typedef struct {
    float x[3];
} tw_fpr;

typedef union f_union {
    float x;
    uint32_t i;
    int32_t ii;
} funion_t;

// void two_prod(float a, float b, float *r, float *e);
void tw_add_sub(const tw_fpr a, const tw_fpr b, tw_fpr *ra, tw_fpr *rb);
void tw_add_sub_ct(const tw_fpr a, const tw_fpr b, tw_fpr *ra, tw_fpr *rb);
tw_fpr tw_sum(const tw_fpr a, const tw_fpr b);
tw_fpr tw_sum_ct(const tw_fpr a, const tw_fpr b);
tw_fpr tw_sum_f_ct(const tw_fpr a, const float b);

tw_fpr tw_prod_acc(const tw_fpr a, const tw_fpr b);
tw_fpr tw_prod_fast(const tw_fpr a, const tw_fpr b);
tw_fpr tw_prod_fast_ct(const tw_fpr a, const tw_fpr b);
tw_fpr tw_prod_fast_f_ct(const float a0, const tw_fpr b);

tw_fpr tw_reci(const tw_fpr x);
tw_fpr tw_reci_fast(const tw_fpr x);
tw_fpr tw_reci_fast_ct(const tw_fpr x);

tw_fpr tw_div(const tw_fpr z, const tw_fpr x);
tw_fpr tw_div_fast(const tw_fpr z, const tw_fpr x);
tw_fpr tw_div_fast_ct(const tw_fpr z, const tw_fpr x);

tw_fpr to_tw(float a, float b, float c);

tw_fpr tw_sqrt(const tw_fpr x);
tw_fpr tw_sqrt_fast(const tw_fpr x);
tw_fpr tw_sqrt_fast_ct(const tw_fpr x);

tw_fpr create_from_u64(uint64_t x);

uint64_t convert_to_u64_round_ct(tw_fpr x);

int32_t tw_trunc_fast(const tw_fpr x);
int64_t tw_trunc_full(const tw_fpr x);
int64_t tw_trunc_full_new(const tw_fpr x);
int32_t tw_rint_fast(const tw_fpr x);
int32_t tw_rint_fast_ct(const tw_fpr x);
int64_t tw_floor_fast(const tw_fpr x);
int32_t tw_floor_fast_ct(const tw_fpr x);
tw_fpr tw_sub(tw_fpr x, tw_fpr y);
tw_fpr tw_sub_ct(tw_fpr x, tw_fpr y);
tw_fpr tw_from_int(int64_t i);
tw_fpr tw_from_int_ct(int32_t i);

#endif
