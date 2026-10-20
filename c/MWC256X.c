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

/* This is MWC256X 1.0, a Marsaglia multiply-with-carry generator with
   period approximately 2^255 whose output is scrambled by X: the result is
   the most recent digit xored with the carry, as in David Thomas's MWC64X
   generator. It is one of the fastest pseudorandom number generators
   available, as its only 128-bit operations are a multiplication and a
   sum, and, since the lag is three, the multiplications of consecutive
   steps are independent: the only dependency between consecutive steps is
   the propagation of the carry (in MWC128X, instead, each multiplication
   depends on the result of the previous step). It is an excellent
   generator based on congruential arithmetic.

   As all MWC generators, the underlying generator simulates a
   multiplicative LCG with prime modulus
   m = 0xff32e7e2bdcb0a32ffffffffffffffffffffffffffffffffffffffffffffffff
   and multiplier given by the inverse of 2^64 modulo m. The modulus has a
   particular form, which forces a defect in the lattice structure: the
   triples of digits x_n, x_{n+1}, x_{n+4} lie on planes at distance about
   2^-32. The defect is unlikely to have any impact in practice: MWC128 has
   a defect of the same size on triples of consecutive digits, which is
   detected only by very large birthday-spacings tests in three dimensions.
   Nonetheless, the scrambler X, which costs a single xor, hides it: such a
   test fails both on the upper and on the lower bits of MWC128, but passes
   on MWC128X, and in our tests on generators of the same kind with 32-bit
   digits, X removed all the lattice structure detectable by
   birthday-spacings tests.

   MWC_A3 maximizes, among about 1.2 million random multipliers in
   [2^64 - 2^56 .. 2^64) for which m is a safe prime, a weighted geometric
   mean of the figures of merit of the spectral test of the underlying
   generator in dimensions 2 to 8.

   For a generator of the same type with stronger theoretical guarantees,
   consider GMWC256X, a Goresky-Klapper generalized multiply-with-carry
   generator whose modulus does not force such a defect.

   MWC256X 1.0 replaces MWC256 1.0: it has a different MWC_A3 and scrambles
   the output with X.
*/

#define MWC_A3 0xff32e7e2bdcb0a33

/* The state must be initialized so that 0 < c < MWC_A3 - 1. If you have a
   64-bit seed, we suggest to set c = 1, and to seed a SplitMix generator
   and use its output to fill x, y, and z: do not use related seeds (e.g.,
   consecutive integers) directly, as they generate correlated sequences.
   To obtain different streams, you can also use the jump functions. */
uint64_t x, y, z, c;

uint64_t next(void) {
	const uint64_t result = z ^ c;
	const __uint128_t t = MWC_A3 * (__uint128_t)x + c;
	x = y;
	y = z;
	z = t;
	c = t >> 64;
	return result;
}


/* The following jump functions use a minimal multiprecision library. They
   read the state as the integer x + y * b + z * b^2 + c * b^3, where
   b = 2^64: since each call to next() multiplies it modulo m by
   A = b^-1 mod m = MWC_A3 * b^2 (as MWC_A3 * b^3 = m + 1), jumping ahead by
   n steps multiplies it by A^n mod m. */

#define MP_SIZE 5
#include "mp.c"

static uint64_t mod[MP_SIZE] = { 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, MWC_A3 - 1 };


/* This is the jump function for the generator. It is equivalent
   to 2^128 calls to next(); it can be used to generate 2^127
   non-overlapping subsequences for parallel computations.

   Equivalent C++ Boost multiprecision code:

   cpp_int b = cpp_int(1) << 64;
   cpp_int b3 = b * b * b;
   cpp_int m = MWC_A3 * b3 - 1;
   cpp_int r = cpp_int("0xa9c5be8f5f3f478b61d0e6496fc072d5d612376996c728e897211409482daa2e");
   cpp_int s = ((x + y * b + z * b * b + c * b3) * r) % m;
   x = uint64_t(s);
   y = uint64_t(s >> 64);
   z = uint64_t(s >> 128);
   c = uint64_t(s >> 192);
*/

void jump(void) {
	static uint64_t jump[MP_SIZE] = { 0x97211409482daa2e, 0xd612376996c728e8, 0x61d0e6496fc072d5, 0xa9c5be8f5f3f478b };
	uint64_t state[MP_SIZE] = { x, y, z, c };
	mp_mul(state, jump, mod);
	x = state[0];
	y = state[1];
	z = state[2];
	c = state[3];
}


/* This is the long-jump function for the generator. It is equivalent to
   2^192 calls to next(); it can be used to generate 2^63 starting points,
   from each of which jump() will generate 2^64 non-overlapping
   subsequences for parallel distributed computations.

   Equivalent C++ Boost multiprecision code:

   cpp_int b = cpp_int(1) << 64;
   cpp_int b3 = b * b * b;
   cpp_int m = MWC_A3 * b3 - 1;
   cpp_int r = cpp_int("0x8ed9eee8dba1c3f30d9fad79f808de021c73ed817e55bf408c9cace7f7c08841");
   cpp_int s = ((x + y * b + z * b * b + c * b3) * r) % m;
   x = uint64_t(s);
   y = uint64_t(s >> 64);
   z = uint64_t(s >> 128);
   c = uint64_t(s >> 192);
*/

void long_jump(void) {
	static uint64_t long_jump[MP_SIZE] = { 0x8c9cace7f7c08841, 0x1c73ed817e55bf40, 0xd9fad79f808de02, 0x8ed9eee8dba1c3f3 };
	uint64_t state[MP_SIZE] = { x, y, z, c };
	mp_mul(state, long_jump, mod);
	x = state[0];
	y = state[1];
	z = state[2];
	c = state[3];
}


/* The following arbitrary-jump functions compute A^n mod m by modular
   exponentiation. */

static uint64_t base[MP_SIZE] = { 0, 0, MWC_A3 };   // A = b^-1 mod m = MWC_A3 * b^2

static void jumpmul(const uint64_t * const r) {
	uint64_t state[MP_SIZE] = { x, y, z, c };
	mp_mul(state, r, mod);
	x = state[0];
	y = state[1];
	z = state[2];
	c = state[3];
}

/* This is the arbitrary-jump function for the generator expressed as c * 2^e.
   It is equivalent to c * 2^e calls to next(); for example, jump_ce(1, 128) is
   equivalent to jump() and jump_ce(1, 192) is equivalent to long_jump().
   Expressing the distance as c * 2^e makes it possible to request both ordinary
   counts (jump_ce(k, 0)) and multiples of power-of-two jumps without handling
   multiple-precision integers.

   Equivalent C++ Boost multiprecision code (carry is the generator's c):

   cpp_int b = cpp_int(1) << 64;
   cpp_int b3 = b * b * b;
   cpp_int m = MWC_A3 * b3 - 1;
   cpp_int A = MWC_A3 * b * b;                   // = b^-1 mod m
   cpp_int r = powm(A, cpp_int(c) << e, m);      // A^(c * 2^e) mod m
   cpp_int s = ((x + y * b + z * b * b + carry * b3) * r) % m;
   x = uint64_t(s);
   y = uint64_t(s >> 64);
   z = uint64_t(s >> 128);
   carry = uint64_t(s >> 192);
*/

void jump_ce(uint64_t c, uint32_t e) {
	uint64_t r[MP_SIZE];
	mp_powmod(r, base, c, mod);
	while (e--) mp_mul(r, r, mod);   // r = (base^c)^(2^e) = base^(c * 2^e) mod m
	jumpmul(r);
}

/* This is the general arbitrary-jump function for the generator. It is
   equivalent to n calls to next(), where n = jump[0] + jump[1] * 2^64 + ... is
   the little-endian integer held in the MP_SIZE - 1 words of jump. Unlike
   jump_ce(), it can express any jump distance. For the jump to be meaningful, n
   should be smaller than the period (approximately 2^255).

   Equivalent C++ Boost multiprecision code (carry is the generator's c):

   cpp_int b = cpp_int(1) << 64;
   cpp_int b3 = b * b * b;
   cpp_int m = MWC_A3 * b3 - 1;
   cpp_int A = MWC_A3 * b * b;                   // = b^-1 mod m
   cpp_int r = powm(A, n, m);                    // A^n mod m
   cpp_int s = ((x + y * b + z * b * b + carry * b3) * r) % m;
   x = uint64_t(s);
   y = uint64_t(s >> 64);
   z = uint64_t(s >> 128);
   carry = uint64_t(s >> 192);
*/

void jump_n(const uint64_t jump[MP_SIZE - 1]) {
	uint64_t r[MP_SIZE];
	mp_powmod_arr(r, base, jump, MP_SIZE - 1, mod);
	jumpmul(r);
}
