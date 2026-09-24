/**
 * @file msl_trig.c
 * The game's own single-precision trigonometry, as the GameCube runs it.
 *
 * sinf and cosf come from MSL (src/MSL/trigf.c with the tables of
 * src/MSL/math_data.c) and atanf from src/melee/lb/lbtrigf.c. The port does
 * not compile trigf.c, and lbtrigf.c's atanf is Metrowerks-only, so these
 * calls used to reach the host CRT, whose results differ from the console's
 * in the last bits on a large share of inputs. atan2f, acosf and asinf
 * (lbtrigf.c, compiled) all go through atanf.
 *
 * Metrowerks compiles the source expressions below to fused PowerPC
 * instructions (fmadds, fnmsubs, fnmadds), which round once. MSVC never
 * fuses, so each of those operations is spelled out with msl_fmaf, an exact
 * single-precision fused multiply-add that needs no FMA hardware. The fused
 * sites follow the retail code: sinf 0x803263D4..0x80326574,
 * cosf 0x80326240..0x803263D0.
 */
#include <stdint.h>
#include <string.h>

/* On 32-bit x86 a function returns its float in the x87 st(0), and MSVC then
 * computes on it there at double precision without rounding to single in
 * between. Every step here has to round exactly as written, so the helpers
 * are inlined and everything stays in SSE registers. */
#define MSL_INLINE static __forceinline

/* Exact single-precision fused multiply-add: a*b + c with one rounding.
 * The product of two floats is exact in double. The sum is formed in double
 * with its rounding error recovered exactly (TwoSum), then rounded to odd:
 * a double with a nonzero error keeps an odd last bit. Rounding that
 * round-to-odd double to float is then correctly rounded, since double has
 * more than two bits beyond float's precision. */
MSL_INLINE float msl_fmaf(float a, float b, float c)
{
    double p = (double) a * (double) b;
    double s = p + (double) c;
    double v = s - p;
    double err = (p - (s - v)) + ((double) c - v);
    uint64_t bits;

    if (err != 0.0 && s - s == 0.0) { /* s - s == 0: s is finite */
        memcpy(&bits, &s, sizeof(bits));
        if ((bits & 1) == 0) {
            /* one ulp toward the discarded error */
            bits += ((err > 0.0) == (s > 0.0)) ? 1 : (uint64_t) -1;
            memcpy(&s, &bits, sizeof(bits));
        }
    }
    return (float) s;
}

/* PowerPC fnmsubs: -(a*b - c), fused. Negating after the subtraction keeps
 * the retail sign of an exact zero. */
MSL_INLINE float msl_fnmsubs(float a, float b, float c)
{
    return -msl_fmaf(a, b, -c);
}

MSL_INLINE uint32_t msl_float_bits(float x)
{
    uint32_t bits;
    memcpy(&bits, &x, sizeof(bits));
    return bits;
}

/* src/MSL/math_data.c */
static const float sincos_on_quadrant[] = { 0, 1, 1, 0, 0, -1, -1, 0 };
static const float sincos_poly[] = {
    0.0000035287617f, 0.0000003089747f, -0.0003259365f,
    -0.00003657235f,  0.015854323f,     0.0024903931f,
    -0.30842513f,     -0.08074551f,     1.0f,
    0.7853982f,
};

/* src/MSL/trigf.c: __four_over_pi_m1, filled by __sinit_trigf_c */
static const float four_over_pi_m1[] = { 0.25f, 0.0232393741608f,
                                         1.70555722434e-7f,
                                         1.86736494323e-11f };

#define MSL_EPSILON 3.45266983e-4f

/* Argument in eighth turns relative to the nearest quadrant, and that
 * quadrant (trigf.c; the four correction terms accumulate fused). */
MSL_INLINE float msl_reduce(float x, int* quadrant)
{
    float z = (2.0f / 3.14159265358979323846f) * x;
    int n = (msl_float_bits(x) & 0x80000000u) ? (int) (z - 0.5f)
                                              : (int) (z + 0.5f);
    float y = msl_fmaf(
        four_over_pi_m1[3], x,
        msl_fmaf(four_over_pi_m1[2], x,
                 msl_fmaf(four_over_pi_m1[1], x,
                          msl_fmaf(four_over_pi_m1[0], x, x - n * 2))));
    *quadrant = n & 3;
    return y;
}

MSL_INLINE float msl_abs(float x)
{
    return x < 0.0f ? -x : x;
}

/* Even polynomial: sincos_poly[0, 2, 4, 6, 8] */
MSL_INLINE float msl_poly_even(float ysq)
{
    return msl_fmaf(
        ysq,
        msl_fmaf(ysq,
                 msl_fmaf(ysq,
                          msl_fmaf(sincos_poly[0], ysq, sincos_poly[2]),
                          sincos_poly[4]),
                 sincos_poly[6]),
        sincos_poly[8]);
}

/* Odd polynomial without the final factor y: sincos_poly[1, 3, 5, 7, 9] */
MSL_INLINE float msl_poly_odd(float ysq)
{
    return msl_fmaf(
        ysq,
        msl_fmaf(ysq,
                 msl_fmaf(ysq,
                          msl_fmaf(sincos_poly[1], ysq, sincos_poly[3]),
                          sincos_poly[5]),
                 sincos_poly[7]),
        sincos_poly[9]);
}

MSL_INLINE float msl_sinf(float x)
{
    int n;
    float y = msl_reduce(x, &n);
    float ysq;

    if (msl_abs(y) < MSL_EPSILON) {
        n <<= 1;
        return msl_fmaf(sincos_poly[9], sincos_on_quadrant[n + 1] * y,
                        sincos_on_quadrant[n]);
    }
    ysq = y * y;
    if (n & 1) {
        n <<= 1;
        return msl_poly_even(ysq) * sincos_on_quadrant[n];
    }
    n <<= 1;
    return (msl_poly_odd(ysq) * y) * sincos_on_quadrant[n + 1];
}

MSL_INLINE float msl_cosf(float x)
{
    int n;
    float y = msl_reduce(x, &n);
    float ysq;

    if (msl_abs(y) < MSL_EPSILON) {
        n <<= 1;
        return msl_fnmsubs(y, sincos_on_quadrant[n], sincos_on_quadrant[n + 1]);
    }
    ysq = y * y;
    if (n & 1) {
        n <<= 1;
        /* retail folds the source negation into the last step (fnmadds) */
        return (-msl_poly_odd(ysq) * y) * sincos_on_quadrant[n];
    }
    n <<= 1;
    return msl_poly_even(ysq) * sincos_on_quadrant[n + 1];
}

float melee_sinf(float x)
{
    return msl_sinf(x);
}

float melee_cosf(float x)
{
    return msl_cosf(x);
}

/* trigf.c: tanf is sin__Ff(x) / cos__Ff(x) */
float melee_tanf(float x)
{
    return msl_sinf(x) / msl_cosf(x);
}

/* src/melee/lb/lbtrigf.c */
static const float atanf_lookup[] = {
    1.0f,
    -0.3333333134651184f,
    0.1999988704919815f,
    -0.14281649887561798f,
    0.11041180044412613f,
    -0.08459755778312683f,
    0.04714243486523628f,
    6.828420162200928f,
    3.239828109741211f,
    2.0f,
    1.4464620351791382f,
    1.1715729236602783f,
    1.039566159248352f,
    7.1350000325764995e-06f,
    8.200000252145401e-07f,
    0.0f,
    6.299999881775875e-07f,
    0.0f,
    0.0f,
    0.0f,
    0.3926900029182434f,
    0.5890486240386963f,
    0.7853981256484985f,
    0.9817469716072083f,
    1.1780970096588135f,
    1.3744460344314575f,
    0.0f,
    9.081698408408556e-06f,
    2.3000000126671694e-08f,
    6.30000016599297e-08f,
    7.040000014058023e-07f,
    2.499999993688107e-07f,
    7.900000014160469e-07f,
    2.414212942123413f,
    1.4966057538986206f,
    1.0f,
    0.6681786179542542f,
    0.4142135679721832f,
    0.1989123672246933f,
    5.620000251838064e-07f,
    0.0f,
    0.0f,
    0.0f,
    0.0f,
    0.0f,
    0.0f,
};

float melee_atanf(float x)
{
    const float silver_ratio = 2.4142136573791504f;
    const float silver_ratio_conjugate = 0.4142135679721832f;
    uint32_t sign = msl_float_bits(x) & 0x80000000u;
    uint32_t bits = msl_float_bits(x) & ~0x80000000u;
    int index = -1;
    int reciprocal = 0;
    float result, squared, poly;

    memcpy(&x, &bits, sizeof(x));
    if (x >= silver_ratio) {
        reciprocal = 1;
        result = 1.0f / x;
    } else if (silver_ratio_conjugate < x) {
        float offset_33, offset_39;
        /* the exponent switch of lbtrigf.c, as thresholds on the bits */
        index = ((int32_t) bits >= 0x3F08D5B9) + ((int32_t) bits >= 0x3F521801) +
                ((int32_t) bits >= 0x3F9BF7EC) + ((int32_t) bits >= 0x3FEF789E);
        offset_39 = atanf_lookup[index + 39];
        offset_33 = atanf_lookup[index + 33];
        result = 1.0f / (offset_33 + (x + offset_39));
        result = msl_fnmsubs(result, atanf_lookup[index + 7], offset_33) +
                 msl_fnmsubs(result, atanf_lookup[index + 13], offset_39);
    } else {
        result = x;
    }

    squared = result * result;
    poly = msl_fmaf(squared, atanf_lookup[6], atanf_lookup[5]);
    poly = msl_fmaf(squared, poly, atanf_lookup[4]);
    poly = msl_fmaf(squared, poly, atanf_lookup[3]);
    poly = msl_fmaf(squared, poly, atanf_lookup[2]);
    poly = msl_fmaf(squared, poly, atanf_lookup[1]);
    result = msl_fmaf(result * squared, poly, result);
    /* index -1 (no reduction) reads the zeros at 26 and 19 */
    result += atanf_lookup[index + 27];
    result += atanf_lookup[index + 20];

    if (reciprocal) {
        result -= 1.57079637f; /* (float) M_PI_2 */
        return sign ? result : -result;
    }
    bits = msl_float_bits(result) | sign;
    memcpy(&result, &bits, sizeof(result));
    return result;
}
