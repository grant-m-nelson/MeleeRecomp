/**
 * @file pc_platform.h
 * Force-included into every translation unit of the PC build
 * (see pc/CMakeLists.txt). Anything the GameCube toolchain provided
 * implicitly, or that MSVC spells differently, is patched up here so the
 * game sources themselves stay untouched wherever possible.
 */
#ifndef PC_PLATFORM_H
#define PC_PLATFORM_H

#define TARGET_PC 1

/* MSVC ABI only (cl.exe or clang-cl). Keep the checks so another toolchain
 * gets a clear error instead of silently compiling something wrong. */
#if !defined(_MSC_VER)
#error "pc_platform.h currently supports MSVC and clang-cl only"
#endif

#if defined(_WIN64) || defined(__LP64__) || (defined(__SIZEOF_POINTER__) && __SIZEOF_POINTER__ != 4)
#error "The PC build must be 32-bit: HSD archives relocate 32-bit pointers in place"
#endif

/* Alignment attributes. The engine asserts 32-byte alignment on buffers it
 * hands to the disc and audio DMA paths, so the attribute has to work.
 * clang-cl accepts the GNU spelling in the trailing position the code uses;
 * MSVC's __declspec(align) cannot go there, so under cl.exe the attribute is
 * compiled out and those asserts fire (build with clang-cl: pc\build.cmd). */
#if defined(__clang__)
#define ATTRIBUTE_ALIGN(num) __attribute__((aligned(num)))
#else
#define ATTRIBUTE_ALIGN(num)
#endif

/* Attribute spellings that <Runtime/platform.h> resolves to GNU syntax. */
#define ATTRIBUTE_NORETURN __declspec(noreturn)
#define DOLPHIN_ATTRIBUTE_NORETURN __declspec(noreturn)
#define UNUSED

/* PowerPC intrinsics used directly by the game/engine sources.
 * __frsqrte is the reciprocal square-root *estimate*; callers refine it with
 * Newton-Raphson steps, so returning the exact value is fine. */
#define __frsqrte(x) (1.0 / sqrt((double) (x)))
#define __fabs(x) fabs((double) (x))
#define __fabsf(x) fabsf(x)

/* MSVC's math.h only exposes M_PI and friends with this defined. */
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#ifndef _CRT_NONSTDC_NO_WARNINGS
#define _CRT_NONSTDC_NO_WARNINGS
#endif

#include <math.h>
#include <stddef.h>
#include <stdint.h>

/* The game ships its own implementations of these (src/melee/lb/lbtrigf.c,
 * lb_00CE.c), matching the GameCube's results. The host CRT defines some of
 * them inline in <math.h>, so rename the game's versions and route every
 * caller to them. */
#define atan2f melee_atan2f
#define acosf melee_acosf
#define asinf melee_asinf
#define expf melee_expf
#define powf melee_powf
float melee_atan2f(float y, float x);
float melee_acosf(float x);
float melee_asinf(float x);
float melee_expf(float x);
float melee_powf(float x, float y);

/* The game's sqrtf (src/MSL/math_ppc.h) is an inline frsqrte + Newton
 * sequence guarded by x > 0: it returns x itself for zero and negative
 * inputs, where the CRT returns NaN. Game code relies on that. The bone
 * dynamics solver (lb_8001044C) takes sqrtf of a negative whenever a link
 * ends within 0.1 of a collider sphere, and with NaN it silently skipped the
 * push away from the body (Fox's tail sinking into a lying fighter, so a hit
 * that connects on the console misses). sqrtf__Ff is the out-of-line MSL
 * function and keeps the CRT's semantics. */
static __inline float pc_msl_sqrtf(float x)
{
    return x > 0.0f ? (sqrtf)(x) : x;
}
#define sqrtf__Ff(x) (sqrtf)(x)
#define sqrtf_accurate(x) pc_msl_sqrtf(x)
#define sqrtf(x) pc_msl_sqrtf(x)

/* sinf, cosf and tanf are MSL's (src/MSL/trigf.c, not compiled here) and
 * atanf is lbtrigf.c's (Metrowerks-only there); pc/src/msl_trig.c has them
 * with the console's fused rounding. Without this they reached the CRT, whose
 * results differ from the console's in the last bits on a large share of
 * inputs; atan2f, acosf and asinf above all call atanf. */
#define sinf melee_sinf
#define cosf melee_cosf
#define tanf melee_tanf
#define atanf melee_atanf
float melee_sinf(float x);
float melee_cosf(float x);
float melee_tanf(float x);
float melee_atanf(float x);

/* Metrowerks setjmp: the game embeds __jmp_buf (a 248-byte PowerPC register
 * image) inside its own structures and calls __setjmp/longjmp on it. On PC we
 * keep the struct shape (so struct layouts and sizes are unchanged) and store
 * the host jmp_buf inside it. See pc/src/pc_setjmp.c. */
#define PC_JMP_BUF_HOST_BYTES 248

#endif /* PC_PLATFORM_H */
