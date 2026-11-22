/*  Written in 2026 by Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide.

Permission to use, copy, modify, and/or distribute this software for any
purpose with or without fee is hereby granted.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR
IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE. */

/*

The step of the generators with 64-bit digits computes a * x + c, where a,
x, and c are 64-bit words, and needs all 128 bits of the result. The code
in this directory computes it using the type __uint128_t, which is not part
of ISO C: GCC and Clang provide it on 64-bit targets, but, for example,
MSVC does not, and nvcc supports 128-bit integers in device code only
from CUDA 11.5 on.

Replacing the 128-bit product with a schoolbook multiplication on 32-bit
halves gives a correct but several times slower step: CPUs compute a
64 x 64 -> 128-bit product with one or two instructions, but compilers do
not recognize the schoolbook multiplication as such. (On GPUs, which
compute 64-bit products using 32-bit instructions anyway, the slowdown is
smaller, but still noticeable.) This header provides a function

	uint64_t mul_add(uint64_t a, uint64_t x, uint64_t c, uint64_t *hi)

returning the lower 64 bits of a * x + c and storing the upper 64 bits in
*hi, which uses, in this order,

- in CUDA device code, unsigned __int128 from CUDA 11.5 on (it is faster
  than __umul64hi() on the GPUs we tested) and __umul64hi() before;
- __uint128_t, if the compiler provides it;
- _umul128() with MSVC on x64 and __umulh() with MSVC on ARM64;
- schoolbook multiplication on 32-bit halves.

For example, the step of MWC128X becomes

	const uint64_t result = x ^ c;
	uint64_t hi;
	x = mul_add(MWC_A1, x, c, &hi);
	c = hi;
	return result;

Defining MUL_ADD_NO_INT128 excludes 128-bit integer types, and defining
MUL_ADD_PORTABLE forces the schoolbook multiplication; both macros are
meant for testing (see test/mul_add.c and test/mul_add.cu).

*/

#ifndef MUL_ADD_H
#define MUL_ADD_H

#include <stdint.h>

#if defined(MUL_ADD_PORTABLE)
#elif defined(__CUDA_ARCH__)
#if !defined(MUL_ADD_NO_INT128) && (__CUDACC_VER_MAJOR__ > 11 || (__CUDACC_VER_MAJOR__ == 11 && __CUDACC_VER_MINOR__ >= 5))
#define MUL_ADD_INT128 unsigned __int128
#else
#define MUL_ADD_UMUL64HI
#endif
#elif defined(__SIZEOF_INT128__) && !defined(MUL_ADD_NO_INT128)
#define MUL_ADD_INT128 __uint128_t
#elif defined(_MSC_VER) && (defined(_M_X64) || defined(_M_AMD64)) && !defined(_M_ARM64EC)
#include <intrin.h>
#pragma intrinsic(_umul128)
#define MUL_ADD_UMUL128
#elif defined(_MSC_VER) && (defined(_M_ARM64) || defined(_M_ARM64EC))
#include <intrin.h>
#define MUL_ADD_UMULH
#endif

#if defined(__CUDACC__)
#define MUL_ADD_INLINE static inline __host__ __device__
#elif defined(_MSC_VER) && !defined(__cplusplus)
#define MUL_ADD_INLINE static __inline
#else
#define MUL_ADD_INLINE static inline
#endif

/* The result fits in 128 bits, as (2^64 - 1)^2 + 2^64 - 1 < 2^128. */
MUL_ADD_INLINE uint64_t mul_add(const uint64_t a, const uint64_t x, const uint64_t c, uint64_t * const hi) {
#if defined(MUL_ADD_INT128)
	const MUL_ADD_INT128 t = (MUL_ADD_INT128)a * x + c;
	*hi = (uint64_t)(t >> 64);
	return (uint64_t)t;
#elif defined(MUL_ADD_UMUL64HI)
	const uint64_t l = a * x + c;
	*hi = __umul64hi(a, x) + (l < c);
	return l;
#elif defined(MUL_ADD_UMUL128)
	uint64_t h;
	const uint64_t l = _umul128(a, x, &h) + c;
	*hi = h + (l < c);
	return l;
#elif defined(MUL_ADD_UMULH)
	const uint64_t l = a * x + c;
	*hi = __umulh(a, x) + (l < c);
	return l;
#else
	/* The halves of c are added to the first two partial products, each of
	   which can absorb two 32-bit addends, as (2^32 - 1)^2 + 2 * (2^32 - 1)
	   = 2^64 - 1. */
	const uint64_t mask = UINT64_C(0xffffffff);
	const uint64_t a0 = a & mask, a1 = a >> 32, x0 = x & mask, x1 = x >> 32;
	const uint64_t t0 = a0 * x0 + (c & mask);
	const uint64_t t1 = a0 * x1 + (t0 >> 32) + (c >> 32);
	const uint64_t t2 = a1 * x0 + (t1 & mask);
	*hi = a1 * x1 + (t1 >> 32) + (t2 >> 32);
	return (t2 << 32) | (t0 & mask);
#endif
}

#endif
