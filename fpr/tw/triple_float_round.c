/*
 * Copyright 2025 NXP
 * SPDX-License-Identifier: MIT
 */

#include "triple_float.h"
#include <fenv.h>
#include <math.h>
#include <stdio.h>

// TODO implement case where x is 0
uint64_t convert_to_u64_round_ct(const tw_fpr x){
    uint64_t r;
    funion_t xt;
    xt.x = x.x[0];
    // set sign
    uint32_t sign = (xt.i >> 31);
    r = (uint64_t)sign << 63;

    uint16_t exponent = ((xt.i >> 23) & 0xFF);

    // mantissa split up in 3 limps, save 31 bits in every limp, total of 93
    // bits of precision (with interleaving zeros) will be rounded to 53 bits of
    // precision in the end
    uint32_t mantissa[3];
    mantissa[0] = 1 << 30;
    mantissa[0] += (xt.i & 0x7FFFFF) << 7;
    mantissa[1] = 0;
    mantissa[2] = 0;

    xt.x = x.x[1];
    uint16_t exponent1 = (xt.i & 0x7FFFFFFF) >> 23;
    uint32_t exponent_diff = exponent - exponent1;
    uint32_t mantissa1_explicit = xt.i & 0x7FFFFF;
    mantissa1_explicit += 1 << 23;

    uint32_t interleaving_zeros = exponent_diff - 24;
    uint32_t sign1 = xt.i >> 31;

    funion_t xt2;
    xt2.x = x.x[2];

    uint16_t exponent2 = (xt2.i & 0x7FFFFFFF) >> 23;
    uint32_t interleaving_zeros2 = exponent1 - exponent2 - 24;
    uint32_t sign2 = xt2.i >> 31;
    uint32_t mantissa2_explicit = xt2.i & 0x7FFFFF;
    mantissa2_explicit += 1 << 23;

    // there are still 7 spaces in mantissa[0],

    int32_t amount_in_0 = 7 - interleaving_zeros;

    uint32_t toggle_add_sub1 = 0u - (sign1 ^ sign) * 2u + 1u;
    uint32_t toggle_add_sub2 = 0u - (sign2 ^ sign) * 2u + 1u;

    uint32_t mask_amount_in_0_smaller_eq_zero =
        ((int32_t)(amount_in_0 - 1) >> 31);
    uint32_t mask_amount_in_0_bigger_zero = ~mask_amount_in_0_smaller_eq_zero;

    uint32_t carry;
    // if (amount_in_0 > 0)
    uint32_t tmp2 = (24 - amount_in_0);

    int32_t left_shift = 14 - interleaving_zeros;

    // Calculate masks for the conditions
    uint32_t mask_left_shift_smaller_eq_zero =
        ((int32_t)(left_shift - 1) >> 31);
    uint32_t mask_left_shift_bigger_zero = ~mask_left_shift_smaller_eq_zero;
    // should only happen if amount_in_0 is smaller or equal to zero
    mask_left_shift_smaller_eq_zero &= mask_amount_in_0_smaller_eq_zero;
    mask_left_shift_bigger_zero &= mask_amount_in_0_smaller_eq_zero;

    // uint32_t tmp_toggle = toggle_add_sub1 & mask_amount_in_0_bigger_zero;

    // operation smaller eq zero
    uint32_t zeros_left_for_2 = interleaving_zeros - 14;
    // mask for tmp < 32
    uint32_t mask1 = (int32_t)(zeros_left_for_2 - 32) >> 31;
    // mask for tmp >= 32 and < 63
    uint32_t mask2 =
        (int32_t)((zeros_left_for_2 - 32) & ~(zeros_left_for_2 - 63)) >> 31;
    mask2 = ~mask2 + 1;

    // should only happen if left_shift is smaller or equal to zero
    mask1 &= mask_left_shift_smaller_eq_zero;
    mask2 &= mask_left_shift_smaller_eq_zero;
    // mask 1
    mantissa[2] += toggle_add_sub1 *
                   (((mantissa1_explicit << (31 - zeros_left_for_2)) &
                     0x7FFFFFFF & mask1) +
                    ((mantissa1_explicit >> (zeros_left_for_2 - 31)) & mask2));

    carry = (mantissa[2] >> 31);

    mantissa[1] +=
        toggle_add_sub1 *
        ((((mantissa1_explicit >> zeros_left_for_2) & mask1) + carry) +
         ((mantissa1_explicit << (31 - tmp2)) & 0x7FFFFFFF &
          mask_amount_in_0_bigger_zero) +
         ((mantissa1_explicit << left_shift) & mask_left_shift_bigger_zero));

    carry = (mantissa[1] >> 31);
    mantissa[0] +=
        toggle_add_sub1 *
        (((mantissa1_explicit >> tmp2) & mask_amount_in_0_bigger_zero) + carry);

    // make sure only 31 bits are used
    mantissa[2] &= 0x7FFFFFFF;
    mantissa[1] &= 0x7FFFFFFF;

    // Calculate spaces_left_in_1 based on the masks
    int32_t spaces_left_in_1 =
        (left_shift & mask_left_shift_bigger_zero) |
        (-zeros_left_for_2 & mask_left_shift_smaller_eq_zero) |
        ((31 - tmp2) & mask_amount_in_0_bigger_zero);
    // endif

    // TODO fix this branching by making it constant time
    uint32_t mask_spaces_left_in_1_smaller_eq_zero =
        ((int32_t)(left_shift - 1) >> 31);
    uint32_t mask_spaces_left_in_1_bigger_zero =
        ~mask_spaces_left_in_1_smaller_eq_zero;

    // if (spaces_left_in_1 > 0)
    uint32_t shift_right1 = 24 - spaces_left_in_1 + interleaving_zeros2;

    // mask for shift_right1 < 32
    int32_t mask_shift_right1_less_32 = (int32_t)(shift_right1 - 32) >> 31;
    mask_shift_right1_less_32 &= mask_spaces_left_in_1_bigger_zero;

    // mask for shift_right1 < 24
    int32_t mask_shift_right1_less_24 = (int32_t)(shift_right1 - 24) >> 31;
    int32_t mask_shift_right1_less_55 = (int32_t)(shift_right1 - 55) >> 31;
    int32_t mask_shift_right1_between_24_inc_55_exc =
        (~mask_shift_right1_less_24) & mask_shift_right1_less_55;

    mask_shift_right1_less_24 &= mask_spaces_left_in_1_bigger_zero;
    mask_shift_right1_between_24_inc_55_exc &=
        mask_spaces_left_in_1_bigger_zero;

    // else
    // 1 is full, some might be in 2, fill 2 with whatever is left
    // shift right is equal to the amount of values filled in 2 and the
    // amount of interleavin zeros
    uint32_t shift_right = (-spaces_left_in_1) + interleaving_zeros2;

    mantissa[2] +=
        toggle_add_sub2 *
        (((mantissa2_explicit >> shift_right) &
          mask_spaces_left_in_1_smaller_eq_zero) +
         (((mantissa2_explicit << (31 - shift_right1)) & 0x7FFFFFFF) &
          mask_shift_right1_less_24) +
         ((mantissa2_explicit >> (shift_right1 - 24)) &
          mask_shift_right1_between_24_inc_55_exc));

    carry = mantissa[2] >> 31;
    mantissa[2] &= 0x7FFFFFFF;

    // // if shift_right1 < 32
    mantissa[1] += toggle_add_sub2 * (((mantissa2_explicit >> shift_right1) &
                                       mask_shift_right1_less_32) +
                                      carry);

    carry = mantissa[1] >> 31;
    mantissa[1] &= 0x7FFFFFFF;
    mantissa[0] += toggle_add_sub2 * carry;
    // endif

    uint32_t exponent_64 = exponent - 127 + 1023;

    // TODO is it possible to get a 1 in bit 32?
    // probably, has to be combined with implicit bit being set because if bit
    // 32 is 1, bit 31 is 0. happens only if all 3 parts of tw_fpr are all 1's,
    // and there are no interleaving zeros, which is incredibly rare. If fixed,
    // test if it works by creating the input

    // if implicit bit is no longer set, the next one will always be set,
    // because not enough will be subtracted to make that one also go to 0
    uint32_t implicit_bit_still_set = mantissa[0] >> 30;
    uint32_t implicit_bit_not_set = implicit_bit_still_set ^ 1;
    exponent_64 -= implicit_bit_not_set;
    mantissa[0] <<= implicit_bit_not_set;
    mantissa[1] <<= implicit_bit_not_set;
    mantissa[0] += mantissa[1] >> 31;
    mantissa[1] &= 0x7FFFFFFF;

    uint32_t rounding = 0;
    uint32_t bit54 = (mantissa[1] & 0b0000000000000000000000100000000) >> 8;
    uint32_t bits_above_54 = (mantissa[1] & 0b11111111) | mantissa[2];
    uint32_t bit53 = (mantissa[1] & 0b0000000000000000000001000000000) >> 9;
    // any_bit_set_above_54 is equal to 1 if any bit in bits_above_54 is set, 0
    // otherwise
    uint32_t any_bit_set_above_54 = bits_above_54 | (bits_above_54 >> 1);
    any_bit_set_above_54 |= any_bit_set_above_54 >> 2;
    any_bit_set_above_54 |= any_bit_set_above_54 >> 4;
    any_bit_set_above_54 |= any_bit_set_above_54 >> 8;
    any_bit_set_above_54 |= any_bit_set_above_54 >> 16;

    rounding = bit54 & (any_bit_set_above_54 | bit53);

    // add rounding to 53th bit
    mantissa[1] += rounding << 9;
    mantissa[0] += mantissa[1] >> 31;
    mantissa[1] &= 0x7FFFFFFF;

    // TODO figure out if it is even possible for this to happen
    // probably not, figure out an input that would trigger this, if it exists
    // test it, if not then this branch can be removed
    if (mantissa[0] >> 30 != 1) {
        exponent_64 -= 1;
        mantissa[0] <<= 1;
        mantissa[1] <<= 1;
    }

    r += (uint64_t)(mantissa[0] & 0x3FFFFFFF) << 22;
    r += mantissa[1] >> 9;
    r += (uint64_t)exponent_64 << 52;

    return r;
}


// input is assumed to be positive, and <= 2^63
int64_t tw_trunc_full(const tw_fpr x) {
    funion_t xt;
    xt.x = x.x[0];

    uint16_t exponent = ((xt.i >> 23) & 0xFF);

    // mantissa split up in 3 limps, save 31 bits in every limp, total of 93
    // bits of precision (with interleaving zeros) will be rounded to 53 bits of
    // precision in the end
    uint32_t mantissa[3];
    mantissa[0] = 1 << 30;
    mantissa[0] += (xt.i & 0x7FFFFF) << 7;
    mantissa[1] = 0;
    mantissa[2] = 0;

    xt.x = x.x[1];
    uint16_t exponent1 = (xt.i & 0x7FFFFFFF) >> 23;
    uint32_t exponent_diff = exponent - exponent1;
    uint32_t mantissa1_explicit = xt.i & 0x7FFFFF;
    mantissa1_explicit += 1 << 23;

    uint32_t interleaving_zeros = exponent_diff - 24;
    uint32_t sign1 = xt.i >> 31;

    funion_t xt2;
    xt2.x = x.x[2];

    uint16_t exponent2 = (xt2.i & 0x7FFFFFFF) >> 23;
    uint32_t interleaving_zeros2 = exponent1 - exponent2 - 24;
    uint32_t sign2 = xt2.i >> 31;
    uint32_t mantissa2_explicit = xt2.i & 0x7FFFFF;
    mantissa2_explicit += 1 << 23;

    // there are still 7 spaces in mantissa[0],

    int32_t amount_in_0 = 7 - interleaving_zeros;

    uint32_t toggle_add_sub1 = - sign1 * 2u + 1u;
    uint32_t toggle_add_sub2 = - sign2 * 2u + 1u;

    uint32_t mask_amount_in_0_smaller_eq_zero =
        ((int32_t)(amount_in_0 - 1) >> 31);
    uint32_t mask_amount_in_0_bigger_zero = ~mask_amount_in_0_smaller_eq_zero;

    uint32_t carry;
    // if (amount_in_0 > 0)
    uint32_t tmp2 = (24 - amount_in_0);

    int32_t left_shift = 14 - interleaving_zeros;

    // Calculate masks for the conditions
    uint32_t mask_left_shift_smaller_eq_zero =
        ((int32_t)(left_shift - 1) >> 31);
    uint32_t mask_left_shift_bigger_zero = ~mask_left_shift_smaller_eq_zero;
    // should only happen if amount_in_0 is smaller or equal to zero
    mask_left_shift_smaller_eq_zero &= mask_amount_in_0_smaller_eq_zero;
    mask_left_shift_bigger_zero &= mask_amount_in_0_smaller_eq_zero;

    // uint32_t tmp_toggle = toggle_add_sub1 & mask_amount_in_0_bigger_zero;

    // operation smaller eq zero
    uint32_t zeros_left_for_2 = interleaving_zeros - 14;
    // mask for tmp < 32
    uint32_t mask1 = (int32_t)(zeros_left_for_2 - 32) >> 31;
    // mask for tmp >= 32 and < 63
    uint32_t mask2 =
        (int32_t)((zeros_left_for_2 - 32) & ~(zeros_left_for_2 - 63)) >> 31;
    mask2 = ~mask2 + 1;

    // should only happen if left_shift is smaller or equal to zero
    mask1 &= mask_left_shift_smaller_eq_zero;
    mask2 &= mask_left_shift_smaller_eq_zero;
    // mask 1
    mantissa[2] += toggle_add_sub1 *
                   (((mantissa1_explicit << (31 - zeros_left_for_2)) &
                     0x7FFFFFFF & mask1) +
                    ((mantissa1_explicit >> (zeros_left_for_2 - 31)) & mask2));

    carry = (mantissa[2] >> 31);

    mantissa[1] +=
        toggle_add_sub1 *
        ((((mantissa1_explicit >> zeros_left_for_2) & mask1) + carry) +
         ((mantissa1_explicit << (31 - tmp2)) & 0x7FFFFFFF &
          mask_amount_in_0_bigger_zero) +
         ((mantissa1_explicit << left_shift) & mask_left_shift_bigger_zero));

    carry = (mantissa[1] >> 31);
    mantissa[0] +=
        toggle_add_sub1 *
        (((mantissa1_explicit >> tmp2) & mask_amount_in_0_bigger_zero) + carry);

    // make sure only 31 bits are used
    mantissa[2] &= 0x7FFFFFFF;
    mantissa[1] &= 0x7FFFFFFF;

    // Calculate spaces_left_in_1 based on the masks
    int32_t spaces_left_in_1 =
        (left_shift & mask_left_shift_bigger_zero) |
        (-zeros_left_for_2 & mask_left_shift_smaller_eq_zero) |
        ((31 - tmp2) & mask_amount_in_0_bigger_zero);
    // endif

    // TODO fix this branching by making it constant time
    uint32_t mask_spaces_left_in_1_smaller_eq_zero =
        ((int32_t)(left_shift - 1) >> 31);
    uint32_t mask_spaces_left_in_1_bigger_zero =
        ~mask_spaces_left_in_1_smaller_eq_zero;

    // if (spaces_left_in_1 > 0)
    uint32_t shift_right1 = 24 - spaces_left_in_1 + interleaving_zeros2;

    // mask for shift_right1 < 32
    int32_t mask_shift_right1_less_32 = (int32_t)(shift_right1 - 32) >> 31;
    mask_shift_right1_less_32 &= mask_spaces_left_in_1_bigger_zero;

    // mask for shift_right1 < 24
    int32_t mask_shift_right1_less_24 = (int32_t)(shift_right1 - 24) >> 31;
    int32_t mask_shift_right1_less_55 = (int32_t)(shift_right1 - 55) >> 31;
    int32_t mask_shift_right1_between_24_inc_55_exc =
        (~mask_shift_right1_less_24) & mask_shift_right1_less_55;

    mask_shift_right1_less_24 &= mask_spaces_left_in_1_bigger_zero;
    mask_shift_right1_between_24_inc_55_exc &=
        mask_spaces_left_in_1_bigger_zero;

    // else
    // 1 is full, some might be in 2, fill 2 with whatever is left
    // shift right is equal to the amount of values filled in 2 and the
    // amount of interleavin zeros
    uint32_t shift_right = (-spaces_left_in_1) + interleaving_zeros2;

    mantissa[2] +=
        toggle_add_sub2 *
        (((mantissa2_explicit >> shift_right) &
          mask_spaces_left_in_1_smaller_eq_zero) +
         (((mantissa2_explicit << (31 - shift_right1)) & 0x7FFFFFFF) &
          mask_shift_right1_less_24) +
         ((mantissa2_explicit >> (shift_right1 - 24)) &
          mask_shift_right1_between_24_inc_55_exc));

    carry = mantissa[2] >> 31;
    mantissa[2] &= 0x7FFFFFFF;

    // // if shift_right1 < 32
    mantissa[1] += toggle_add_sub2 * (((mantissa2_explicit >> shift_right1) &
                                       mask_shift_right1_less_32) +
                                      carry);

    carry = mantissa[1] >> 31;
    mantissa[1] &= 0x7FFFFFFF;
    mantissa[0] += toggle_add_sub2 * carry;
    // endif

    // uint32_t exponent_64 = exponent - 127 + 1023;

    // TODO is it possible to get a 1 in bit 32?
    // probably, has to be combined with implicit bit being set because if bit
    // 32 is 1, bit 31 is 0. happens only if all 3 parts of tw_fpr are all 1's,
    // and there are no interleaving zeros, which is incredibly rare. If fixed,
    // test if it works by creating the input

    // if implicit bit is no longer set, the next one will always be set,
    // because not enough will be subtracted to make that one also go to 0
    uint32_t implicit_bit_still_set = mantissa[0] >> 30;
    uint32_t implicit_bit_not_set = implicit_bit_still_set ^ 1;
    exponent -= implicit_bit_not_set;
    mantissa[0] <<= implicit_bit_not_set;
    mantissa[1] <<= implicit_bit_not_set;
    mantissa[0] += mantissa[1] >> 31;
    mantissa[1] &= 0x7FFFFFFF;

    // we assume input is always positive, and lower then 62
    uint32_t e = (127 + 62) - exponent;

    uint64_t r = (uint64_t)mantissa[0] << 32;
    r += (uint64_t)mantissa[1] << 1;
    r += (uint64_t)mantissa[2] >> 30;
    r >>= e;


    uint64_t mask = ((int64_t)(exponent - 127) >> 31);// | ((int64_t)(253 - exponent) >> 31);

    return ~mask & r;
}
// not working
int64_t tw_trunc_full_new(const tw_fpr x) {
    fesetround(FE_TOWARDZERO);
    float xf = x.x[0] + x.x[1];
    uint64_t xi = xf;

    if (xf <= 0x1e24) {

        xf = (uint64_t)(x.x[1] + x.x[2]);
        fesetround(FE_TONEAREST);
        // if(xf > 0)
        xi += (uint64_t)xf;
    }

    fesetround(FE_TONEAREST);
    // uint64_t xi = xf;
    return xi;
    // return x.x[0];
}

/**
 * Only works if inputed floating point value is in the range [2^-24, 2^23]
 */
int32_t tw_trunc_fast(const tw_fpr x) {
    fesetround(FE_TOWARDZERO);
    float xf = x.x[0] + x.x[1];
    fesetround(FE_TONEAREST);

    uint32_t xi = xf;
    return xi;
}

/**
 * Only works if inputed floating point value is in the range (2^-23, 2^23)
 * still has to be made constant time
 */

int32_t tw_rint_fast(const tw_fpr x) {
    float fractional_part = x.x[0] - floorf(x.x[0]);

    int32_t xi = 0;
    if (fractional_part == 0.5f) {
        xi = x.x[1] > 0 ? (int)ceilf(x.x[0]) : (int)floorf(x.x[0]);
    } else {
        xi = rintf(x.x[0]);
    }

    // xu.x = x.x[1];
    return xi;
}

// Extract mantissa and normalize it to 24 bits (23 bits + implicit 1)
static inline int32_t float_rint(const uint32_t x) {

    uint32_t m = ((x << 7) | ((uint32_t)1 << 30)) & (((uint32_t)1 << 31) - 1);
    int32_t e = (127 + 30) - ((int32_t)(x >> 23) & 0xFF);

    // Handle very small values (e.g., denormals or zero)
    m &= (uint32_t)((int32_t)(e - 32) >> 16);
    e &= 31;

    // Prepare for rounding
    uint32_t z = m << (31 - e);
    uint32_t y = ((z & 0x3FFFFFFF) + 0x3FFFFFFF) >> 1;
    uint32_t cc = (0xC8 >> (unsigned)((z | y) >> 29)) & 1;

    // Apply rounding
    m = (m >> e) + cc;

    // Apply sign
    uint32_t s = (uint32_t)(*(int32_t *)&x >> 31);
    m = (m ^ s) - s;

    return (int32_t)m;
}
static inline int32_t float_floor(const uint32_t x){

    /* We extract the mantissa as in fpr_rint(), but then we apply the
       sign bit to it; truncation from the integer shift will then
       yield the proper result. */
    uint32_t m = ((x << 7) | ((uint32_t)1 << 30)) & (((uint32_t)1 << 31) - 1);
    uint32_t s = (uint32_t)(*(int32_t *)&x >> 31);
    m = (m ^ s) - s;

    /* Get the shift count. */
    int32_t e = (127 + 30) - ((int32_t)(x >> 23) & 0xFF);

    /* If the shift count is 64 or more, then the value should be 0
       or -1, depending on the sign bit. Note that if the source is
       "minus zero" then we round to -1. We only need to saturate the
       shift count at 63. */
    uint32_t ue = (uint32_t)e;
    ue = (ue | ((31 - ue) >> 16)) & 31;
    // ue = (ue | ((31 - ue) >> 16)) & 31;
    // return 0;
    return ((int32_t)m) >> ue;
}

/**
 * Only works if inputed floating point value is in the range (2^-23, 2^23)
 * still has to be made constant time
 */
int32_t tw_rint_fast_ct(const tw_fpr x) {
    funion_t xff;
    xff.x = x.x[0];
    float fractional_part = x.x[0] - float_floor(xff.i);
    funion_t xf;
    xf.x = fractional_part;

    uint32_t zero_point_5 = 0x3F000000;

    // all 0's if eq, 1 if neq
    int mask_halfway_point =
        (int32_t)((xf.i - zero_point_5) | (zero_point_5 - xf.i)) >> 31;

    funion_t xf1;
    xf1.x = x.x[1];

    int32_t round, trunc;

    // round = rint_f32(xff.i);
    round = float_rint(xff.i);
    trunc = __builtin_truncf(x.x[0]);

    int32_t xi = (round & mask_halfway_point) | (trunc & ~mask_halfway_point);
    int sign1 = xff.i >> 31;
    int sign2 = xf1.i >> 31;
    xi += (1 - 2 * (sign1 | sign2)) & ~mask_halfway_point;

    return xi;
}

/**
 * Only works if inputed floating point value is in the range [2^-24, 2^23]
 * still has to be made constant time
 */
int64_t tw_floor_fast(const tw_fpr x) {
    fesetround(FE_DOWNWARD);

    float xf = x.x[0] + x.x[1];

    fesetround(FE_TONEAREST);

    float xi = floorf(xf);
    return (int64_t)xi;
}


int32_t tw_floor_fast_ct(const tw_fpr x){
    funion_t xf;
    fesetround(FE_DOWNWARD);

    xf.x = x.x[0] + x.x[1];

    fesetround(FE_TONEAREST);

    return float_floor(xf.i);
}

tw_fpr tw_from_int_ct(int32_t i) {
    float x0 = (float)i;
    return (tw_fpr){{x0, 0, 0}};
}

tw_fpr tw_from_int(int64_t i) {
    if (i == 0) {
        return (tw_fpr){{0, 0, 0}};
    }
    if (i < 16777216 && i > -16777216) {
        float x0 = (float)i;
        return (tw_fpr){{x0, 0, 0}};
    } else {
        printf("not possible yet\n");
        return (tw_fpr){{0, 0, 0}};
    }
}
