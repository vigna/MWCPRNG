/*  Written in 2021-2026 by Sebastiano Vigna (vigna@acm.org)

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

#include <stdint.h>
#include <string.h>

/* This is GMWC128X 1.0, a Goresky-Klapper generalized multiply-with-carry
   generator (see their paper "Efficient Multiply-with-Carry Random Number
   Generators with Maximal Period", ACM Trans. Model. Comput. Simul., 13(4), pp.
   310-321, 2003) with period approximately 2^127 whose output is scrambled by
   X, a generalization of the scrambler of David Thomas's MWC64X generator: the
   result is the next digit xored with the upper 64 bits of the 128-bit value
   t = GMWC_A1 * x + c from which the next digit and carry are computed. In an
   MWC generator these bits are the next carry, so X degenerates to Thomas's
   scrambler, the digit xored with the carry; in a GMWC generator the digit
   xored with the carry would be biased, as the range of the carry depends on
   the digit. While in general slower than a scrambled linear generator, it is
   an excellent generator based on congruential arithmetic.

   As all MWC generators, the underlying generator simulates a multiplicative
   LCG with prime modulus m = 0xfff0cfb4ed6366fa0000000068935f57 and multiplier
   given by the inverse of 2^64 modulo m. The major difference with Marsaglia's
   MWC generators is that the modulus has a more general form: this additional
   freedom avoids the defect in the lattice structure that the particular form
   of the modulus forces on MWC generators (see MWC128X), at the price of an
   additional 64-bit multiplication, an additional 128-bit multiplication, and
   an additional 128-bit sum. Both additional multiplications lie on the
   dependency chain of the carry, independently of the lag, so GMWC generators
   are slower than MWC generators. The scrambler X costs a single xor: in our
   tests on generators of the same kind with 32-bit digits, it removed all the
   lattice structure detectable by birthday-spacings tests.

   GMWC_A1 and GMWC_MINUSA0 maximize, among about 11 million random candidates
   with GMWC_A1 in [2^64 - 2^56 .. 2^64) and odd GMWC_MINUSA0 < 2^32 for which m
   is a safe prime, a weighted geometric mean of the figures of merit of the
   spectral test of the underlying generator in dimensions 2 to 8.

   GMWC128X 1.0 replaces GMWC128 1.0: it has different GMWC_A1 and GMWC_MINUSA0,
   and scrambles the output with X.
*/

#define GMWC_MINUSA0 0x68935f57
#define GMWC_A0INV 0x8316d3ad5de20399
#define GMWC_A1 0xfff0cfb4ed6366fa

/* The state must be initialized so that GMWC_MINUSA0 <= c <= GMWC_A1 (a
   sufficient condition: valid states may have carries outside this range).
   If you have a 64-bit seed, we suggest to set c = 2 * GMWC_MINUSA0 (with
   c = GMWC_MINUSA0 or c = GMWC_A1 the states obtained from different seeds
   are arithmetically related), and to seed a SplitMix generator and use
   its output to fill x: do not use related seeds (e.g., consecutive
   integers) directly, as they generate correlated sequences. To obtain
   different streams, you can also use the jump functions. */
uint64_t x, c;

uint64_t next(void) {
	const __uint128_t t = GMWC_A1 * (__uint128_t)x + c;
	x = GMWC_A0INV * (uint64_t)t;
	c = (t + GMWC_MINUSA0 * (__uint128_t)x) >> 64;
	return x ^ (uint64_t)(t >> 64);
}


/* The following jump functions use a minimal multiprecision library. They
   map the state to the integer h = c * b - GMWC_MINUSA0 * x mod m, where
   b = 2^64: since each call to next() multiplies h modulo m by
   A = b^-1 mod m, jumping ahead by n steps multiplies h by A^n mod m. The
   state is then recovered as x = GMWC_A0INV * h mod b, where
   GMWC_A0INV = (-GMWC_MINUSA0)^-1 mod b, and
   c = (h + GMWC_MINUSA0 * x) / b. */

#define MP_SIZE 3
#include "mp.c"

static uint64_t mod[MP_SIZE] = { GMWC_MINUSA0, GMWC_A1 };
static __uint128_t mod128 = ((__uint128_t)GMWC_A1 << 64) + GMWC_MINUSA0;


/* This is the jump function for the generator. It is equivalent
   to 2^64 calls to next(); it can be used to generate 2^63
   non-overlapping subsequences for parallel computations.

   Equivalent C++ Boost multiprecision code:

   cpp_int b = cpp_int(1) << 64;
   cpp_int m = GMWC_A1 * b + GMWC_MINUSA0;
   cpp_int r = cpp_int("0x60e680f9fb32a4e53d3454cf128d6951");
   cpp_int t = ((m + c * b - GMWC_MINUSA0 * cpp_int(x)) * r) % m;
   x = uint64_t(t * GMWC_A0INV);
   c = uint64_t((t + GMWC_MINUSA0 * cpp_int(x)) >> 64);
*/

void jump(void) {
	static uint64_t jump[MP_SIZE] = { 0x3d3454cf128d6951, 0x60e680f9fb32a4e5 };
	const __uint128_t cc = (__uint128_t)c << 64;
	const __uint128_t xx = GMWC_MINUSA0 * (__uint128_t)x;
	const __uint128_t s = cc >= xx ? cc - xx : (mod128 - xx) + cc;
	uint64_t state[MP_SIZE] = { (uint64_t)s, (uint64_t)(s >> 64) };
	mp_mul(state, jump, mod);
	x = GMWC_A0INV * state[0];
	c = ((((__uint128_t)state[1] << 64) + state[0]) + GMWC_MINUSA0 * (__uint128_t)x) >> 64;
}


/* This is the long-jump function for the generator. It is equivalent to
   2^96 calls to next(); it can be used to generate 2^31 starting points,
   from each of which jump() will generate 2^32 non-overlapping
   subsequences for parallel distributed computations.

   Equivalent C++ Boost multiprecision code:

   cpp_int b = cpp_int(1) << 64;
   cpp_int m = GMWC_A1 * b + GMWC_MINUSA0;
   cpp_int r = cpp_int("0xefbb70ae908b3b0e3d469e7c688d3f2f");
   cpp_int t = ((m + c * b - GMWC_MINUSA0 * cpp_int(x)) * r) % m;
   x = uint64_t(t * GMWC_A0INV);
   c = uint64_t((t + GMWC_MINUSA0 * cpp_int(x)) >> 64);
*/

void long_jump(void) {
	static uint64_t long_jump[MP_SIZE] = { 0x3d469e7c688d3f2f, 0xefbb70ae908b3b0e };
	const __uint128_t cc = (__uint128_t)c << 64;
	const __uint128_t xx = GMWC_MINUSA0 * (__uint128_t)x;
	const __uint128_t s = cc >= xx ? cc - xx : (mod128 - xx) + cc;
	uint64_t state[MP_SIZE] = { (uint64_t)s, (uint64_t)(s >> 64) };
	mp_mul(state, long_jump, mod);
	x = GMWC_A0INV * state[0];
	c = ((((__uint128_t)state[1] << 64) + state[0]) + (__uint128_t)GMWC_MINUSA0 * x) >> 64;
}


/* The following arbitrary-jump functions compute A^n mod m by modular
   exponentiation. Since m has no special form, A is given explicitly. */

static uint64_t base[MP_SIZE] = { 0x02da9505fade3656, 0x830f0c9c3dd132cf };   // A = b^-1 mod m

static void jumpmul(const uint64_t * const r) {
	const __uint128_t cc = (__uint128_t)c << 64;
	const __uint128_t xx = GMWC_MINUSA0 * (__uint128_t)x;
	const __uint128_t s = cc >= xx ? cc - xx : (mod128 - xx) + cc;
	uint64_t state[MP_SIZE] = { (uint64_t)s, (uint64_t)(s >> 64) };
	mp_mul(state, r, mod);
	x = GMWC_A0INV * state[0];
	c = ((((__uint128_t)state[1] << 64) + state[0]) + GMWC_MINUSA0 * (__uint128_t)x) >> 64;
}

/* This is the arbitrary-jump function for the generator expressed as c * 2^e.
   It is equivalent to c * 2^e calls to next(); for example, jump_ce(1, 64) is
   equivalent to jump() and jump_ce(1, 96) is equivalent to long_jump().
   Expressing the distance as c * 2^e makes it possible to request both ordinary
   counts (jump_ce(k, 0)) and multiples of power-of-two jumps without handling
   multiple-precision integers.

   Equivalent C++ Boost multiprecision code (carry is the generator's c):

   cpp_int b = cpp_int(1) << 64;
   cpp_int m = GMWC_A1 * b + GMWC_MINUSA0;
   cpp_int A("0x830f0c9c3dd132cf02da9505fade3656");   // = b^-1 mod m
   cpp_int r = powm(A, cpp_int(c) << e, m);           // A^(c * 2^e) mod m
   cpp_int t = ((m + carry * b - GMWC_MINUSA0 * cpp_int(x)) * r) % m;
   x = uint64_t(t * GMWC_A0INV);
   carry = uint64_t((t + GMWC_MINUSA0 * cpp_int(x)) >> 64);
*/

void jump_ce(uint64_t c, uint32_t e) {
	uint64_t r[MP_SIZE];
	mp_powmod(r, base, c, mod);
	while (e--) mp_mul(r, r, mod);   // r = (base^c)^(2^e) = base^(c * 2^e) mod m
	jumpmul(r);
}

/* This is the general arbitrary-jump function for the generator. It is
   equivalent to n calls to next(), where n = jump[0] + jump[1] * 2^64 is the
   little-endian integer held in the MP_SIZE - 1 words of jump. Unlike
   jump_ce(), it can express any jump distance. For the jump to be meaningful, n
   should be smaller than the period (approximately 2^127).

   Equivalent C++ Boost multiprecision code (carry is the generator's c):

   cpp_int b = cpp_int(1) << 64;
   cpp_int m = GMWC_A1 * b + GMWC_MINUSA0;
   cpp_int A("0x830f0c9c3dd132cf02da9505fade3656");   // = b^-1 mod m
   cpp_int r = powm(A, n, m);                         // A^n mod m
   cpp_int t = ((m + carry * b - GMWC_MINUSA0 * cpp_int(x)) * r) % m;
   x = uint64_t(t * GMWC_A0INV);
   carry = uint64_t((t + GMWC_MINUSA0 * cpp_int(x)) >> 64);
*/

void jump_n(const uint64_t jump[MP_SIZE - 1]) {
	uint64_t r[MP_SIZE];
	mp_powmod_arr(r, base, jump, MP_SIZE - 1, mod);
	jumpmul(r);
}
