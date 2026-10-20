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

#include <stdint.h>
#include <string.h>

/* This is MWC128X/32 1.0, a Marsaglia multiply-with-carry generator with 32-bit
   output and period approximately 2^127 whose output is scrambled by X: the
   result is the most recent digit xored with the carry, as in David Thomas's
   MWC64X generator. Its only 64-bit operations are a multiplication and a sum,
   and, since the lag is three, the multiplications of consecutive steps are
   independent: the only dependency between consecutive steps is the propagation
   of the carry (in MWC64X/32, instead, each multiplication depends on the
   result of the previous step). It is our suggested generator for 32-bit
   outputs.

   As all MWC generators, the underlying generator simulates a multiplicative
   LCG with prime modulus m = 0xffcccc65ffffffffffffffffffffffff and multiplier
   given by the inverse of 2^32 modulo m. The modulus has a particular form,
   which forces a defect in the lattice structure: the triples of digits x_n,
   x_{n+1}, x_{n+4} lie on planes at distance about 2^-16, and a
   three-dimensional birthday-spacings test fails in a few seconds both on the
   upper and on the lower bits. The scrambler X, which costs a single xor, hides
   it: in our tests, X removed all the lattice structure detectable by
   birthday-spacings tests. Do not use the unscrambled digits.

   MWC_A3 maximizes, among all 5601 multipliers in [2^32 - 2^24 .. 2^32) for
   which m is a safe prime, a weighted geometric mean of the figures of merit of
   the spectral test of the underlying generator in dimensions 2 to 8.
*/

#define MWC_A3 0xffcccc66

/* The state must be initialized so that 0 < c < MWC_A3 - 1. If you have a
   64-bit seed, we suggest to set c = 1, and to seed a SplitMix generator
   and use its output to fill x, y, and z: do not use related seeds (e.g.,
   consecutive integers) directly, as they generate correlated sequences.
   To obtain different streams, you can also use the jump functions. */
uint32_t x, y, z, c;

uint32_t next(void) {
	const uint32_t result = z ^ c;
	const uint64_t t = MWC_A3 * (uint64_t)x + c;
	x = y;
	y = z;
	z = t;
	c = t >> 32;
	return result;
}


/* The following jump functions use a minimal multiprecision library. They
   read the state as the integer x + y * b + z * b^2 + c * b^3, where
   b = 2^32, stored in two 64-bit words: since each call to next()
   multiplies it modulo m by A = b^-1 mod m = MWC_A3 * b^2 (as
   MWC_A3 * b^3 = m + 1), jumping ahead by n steps multiplies it by
   A^n mod m. */

#define MP_SIZE 3
#include "mp.c"

static uint64_t mod[MP_SIZE] = { 0xffffffffffffffff, (uint64_t)(MWC_A3 - 1) << 32 | 0xffffffff };

static void jumpmul(const uint64_t * const r) {
	uint64_t state[MP_SIZE] = { (uint64_t)y << 32 | x, (uint64_t)c << 32 | z };
	mp_mul(state, r, mod);
	x = state[0];
	y = state[0] >> 32;
	z = state[1];
	c = state[1] >> 32;
}


/* This is the jump function for the generator. It is equivalent
   to 2^64 calls to next(); it can be used to generate 2^63
   non-overlapping subsequences for parallel computations.

   Equivalent C++ Boost multiprecision code:

   cpp_int b = cpp_int(1) << 32;
   cpp_int b3 = b * b * b;
   cpp_int m = MWC_A3 * b3 - 1;
   cpp_int r = cpp_int("0xe67750f4f6c4e62f3df51ad941ad5cca");
   cpp_int s = ((x + y * b + z * b * b + c * b3) * r) % m;
   x = uint32_t(s);
   y = uint32_t(s >> 32);
   z = uint32_t(s >> 64);
   c = uint32_t(s >> 96);
*/

void jump(void) {
	static uint64_t jump[MP_SIZE] = { 0x3df51ad941ad5cca, 0xe67750f4f6c4e62f };
	jumpmul(jump);
}


/* This is the long-jump function for the generator. It is equivalent to
   2^96 calls to next(); it can be used to generate 2^31 starting points,
   from each of which jump() will generate 2^32 non-overlapping
   subsequences for parallel distributed computations.

   Equivalent C++ Boost multiprecision code:

   cpp_int b = cpp_int(1) << 32;
   cpp_int b3 = b * b * b;
   cpp_int m = MWC_A3 * b3 - 1;
   cpp_int r = cpp_int("0xfef31f136713bbb48cfada360d80d5a5");
   cpp_int s = ((x + y * b + z * b * b + c * b3) * r) % m;
   x = uint32_t(s);
   y = uint32_t(s >> 32);
   z = uint32_t(s >> 64);
   c = uint32_t(s >> 96);
*/

void long_jump(void) {
	static uint64_t long_jump[MP_SIZE] = { 0x8cfada360d80d5a5, 0xfef31f136713bbb4 };
	jumpmul(long_jump);
}


/* The following arbitrary-jump functions compute A^n mod m by modular
   exponentiation. */

static uint64_t base[MP_SIZE] = { 0, MWC_A3 };   // A = b^-1 mod m = MWC_A3 * b^2 = MWC_A3 * 2^64

/* This is the arbitrary-jump function for the generator expressed as c * 2^e.
   It is equivalent to c * 2^e calls to next(); for example, jump_ce(1, 64) is
   equivalent to jump() and jump_ce(1, 96) is equivalent to long_jump().
   Expressing the distance as c * 2^e makes it possible to request both ordinary
   counts (jump_ce(k, 0)) and multiples of power-of-two jumps without handling
   multiple-precision integers.

   Equivalent C++ Boost multiprecision code (carry is the generator's c):

   cpp_int b = cpp_int(1) << 32;
   cpp_int b3 = b * b * b;
   cpp_int m = MWC_A3 * b3 - 1;
   cpp_int A = MWC_A3 * b * b;                   // = b^-1 mod m
   cpp_int r = powm(A, cpp_int(c) << e, m);      // A^(c * 2^e) mod m
   cpp_int s = ((x + y * b + z * b * b + carry * b3) * r) % m;
   x = uint32_t(s);
   y = uint32_t(s >> 32);
   z = uint32_t(s >> 64);
   carry = uint32_t(s >> 96);
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

   cpp_int b = cpp_int(1) << 32;
   cpp_int b3 = b * b * b;
   cpp_int m = MWC_A3 * b3 - 1;
   cpp_int A = MWC_A3 * b * b;                   // = b^-1 mod m
   cpp_int r = powm(A, n, m);                    // A^n mod m
   cpp_int s = ((x + y * b + z * b * b + carry * b3) * r) % m;
   x = uint32_t(s);
   y = uint32_t(s >> 32);
   z = uint32_t(s >> 64);
   carry = uint32_t(s >> 96);
*/

void jump_n(const uint64_t jump[MP_SIZE - 1]) {
	uint64_t r[MP_SIZE];
	mp_powmod_arr(r, base, jump, MP_SIZE - 1, mod);
	jumpmul(r);
}
