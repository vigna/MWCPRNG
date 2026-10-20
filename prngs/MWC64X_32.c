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

/* This is MWC64X/32 1.0, a Marsaglia multiply-with-carry generator with 32-bit
   output and period approximately 2^63 whose output is scrambled by X: the
   result is the digit xored with the carry. It is David Thomas's MWC64X
   generator with a different multiplier. Its only 64-bit operations are a
   multiplication and a sum, but, since the lag is one, each multiplication
   depends on the result of the previous step (in MWC128X/32, instead, the
   multiplications of consecutive steps are independent). It is meant for
   environments in which registers are scarce, such as GPUs: otherwise, use
   MWC128X/32, which has a larger period and is at least as fast.

   As all MWC generators, the underlying generator simulates a multiplicative
   LCG with prime modulus m = 0xfff145f7ffffffff and multiplier given by the
   inverse of 2^32 modulo m. The modulus has a particular form, which forces a
   defect in the lattice structure: triples of consecutive digits lie on planes
   at distance about 2^-16, and a three-dimensional birthday-spacings test fails
   in a few seconds both on the upper and on the lower bits. Moreover, as for
   every congruential generator with a 64-bit modulus, the lattice structure of
   some other projections can be detected quickly. The scrambler X, which costs
   a single xor, hides all of this: in our tests, X removed all the lattice
   structure detectable by birthday-spacings tests. Do not use the unscrambled
   digits.

   MWC_A1 maximizes, among all 22982 multipliers in [2^32 - 2^24 .. 2^32) for
   which m is a safe prime, a weighted geometric mean of the figures of merit of
   the spectral test of the underlying generator in dimensions 2 to 8.
*/

#define MWC_A1 0xfff145f8

/* The state must be initialized so that 0 < c < MWC_A1 - 1. If you have a
   64-bit seed, we suggest to set c = 1, and to seed a SplitMix generator
   and use its output to fill x: do not use related seeds (e.g.,
   consecutive integers) directly, as they generate correlated sequences.
   To obtain different streams, you can also use the jump functions. */
uint32_t x, c;

uint32_t next(void) {
	const uint32_t result = x ^ c;
	const uint64_t t = MWC_A1 * (uint64_t)x + c;
	x = t;
	c = t >> 32;
	return result;
}


/* The modulus m = MWC_A1 * 2^32 - 1 fits in 64 bits, so the following jump
   functions need only 128-bit arithmetic. They read the state as the
   integer x + c * b, where b = 2^32: since each call to next() multiplies
   it modulo m by A = b^-1 mod m = MWC_A1 (as MWC_A1 * b = m + 1), jumping
   ahead by n steps multiplies it by A^n mod m. */

static const uint64_t mod = ((uint64_t)MWC_A1 << 32) - 1;

static uint64_t mulmod(const uint64_t a, const uint64_t b) {
	return (__uint128_t)a * b % mod;
}

static void jumpmul(const uint64_t r) {
	const uint64_t s = mulmod((uint64_t)c << 32 | x, r);
	x = s;
	c = s >> 32;
}


/* This is the jump function for the generator. It is equivalent
   to 2^32 calls to next(); it can be used to generate 2^31
   non-overlapping subsequences for parallel computations. */

void jump(void) {
	jumpmul(0x1012f7168d155184);
}


/* This is the long-jump function for the generator. It is equivalent to
   2^48 calls to next(); it can be used to generate 2^15 starting points,
   from each of which jump() will generate 2^16 non-overlapping
   subsequences for parallel distributed computations. */

void long_jump(void) {
	jumpmul(0x7b046c7221442678);
}


/* The following arbitrary-jump functions compute A^n mod m by modular
   exponentiation. */

static uint64_t powmod(uint64_t n) {
	uint64_t r = 1, a = MWC_A1;
	for (; n != 0; n >>= 1, a = mulmod(a, a)) if (n & 1) r = mulmod(r, a);
	return r;
}

/* This is the arbitrary-jump function for the generator expressed as c * 2^e.
   It is equivalent to c * 2^e calls to next(); for example, jump_ce(1, 32) is
   equivalent to jump() and jump_ce(1, 48) is equivalent to long_jump().
   Expressing the distance as c * 2^e makes it possible to request both ordinary
   counts (jump_ce(k, 0)) and multiples of power-of-two jumps. */

void jump_ce(uint64_t c, uint32_t e) {
	uint64_t r = powmod(c);
	while (e--) r = mulmod(r, r);   // r = (A^c)^(2^e) = A^(c * 2^e) mod m
	jumpmul(r);
}

/* This is the general arbitrary-jump function for the generator. It is
   equivalent to n calls to next(). Unlike jump_ce(), it can express any jump
   distance. For the jump to be meaningful, n should be smaller than the period
   (approximately 2^63). */

void jump_n(const uint64_t n) {
	jumpmul(powmod(n));
}
