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

/* This is GMWC256X 1.0, a Goresky-Klapper generalized multiply-with-carry
   generator (see their paper "Efficient Multiply-with-Carry Random Number
   Generators with Maximal Period", ACM Trans. Model. Comput. Simul., 13(4), pp.
   310-321, 2003) with period approximately 2^255 whose output is scrambled by
   X, a generalization of the scrambler of David Thomas's MWC64X generator: the
   result is the next digit xored with the upper 64 bits of the 128-bit value
   t = GMWC_A3 * x + c from which the next digit and carry are computed. In an
   MWC generator these bits are the next carry, so X degenerates to Thomas's
   scrambler, the digit xored with the carry; in a GMWC generator the digit
   xored with the carry would be biased, as the range of the carry depends on
   the digit. While in general slower than a scrambled linear generator, it is
   an excellent generator based on congruential arithmetic.

   As all MWC generators, the underlying generator simulates a multiplicative
   LCG with prime modulus
   m = 0xff76e600f25ae046000000000000000000000000000000000034dabaf938e45b and
   multiplier given by the inverse of 2^64 modulo m. The major difference with
   Marsaglia's MWC generators is that the modulus has a more general form: this
   additional freedom weakens considerably the defect in the lattice structure
   that the particular form of the modulus forces on MWC generators (see
   MWC256X), at the price of an additional 64-bit multiplication, an additional
   128-bit multiplication, and an additional 128-bit sum. Both additional
   multiplications lie on the dependency chain of the carry, independently of
   the lag, so GMWC generators are slower than MWC generators. The scrambler X
   costs a single xor: in our tests on generators of the same kind with 32-bit
   digits, it removed all the lattice structure detectable by birthday-spacings
   tests.

   GMWC_A3 and GMWC_MINUSA0 maximize, among about 1.3 million random candidates
   with GMWC_A3 in [2^64 - 2^56 .. 2^64) and odd GMWC_MINUSA0 < 2^56 for which m
   is a safe prime, a weighted geometric mean of the figures of merit of the
   spectral test of the underlying generator in dimensions 2 to 8.

   GMWC256X 1.0 replaces GMWC256 1.0: it has different GMWC_A3 and GMWC_MINUSA0,
   and scrambles the output with X.
*/

#define GMWC_MINUSA0 0x34dabaf938e45b
#define GMWC_A0INV 0xac8f71df86d6542d
#define GMWC_A3 0xff76e600f25ae046

/* The state must be initialized so that GMWC_MINUSA0 <= c <= GMWC_A3 (a
   sufficient condition: valid states may have carries outside this range).
   If you have a 64-bit seed, we suggest to set c = 2 * GMWC_MINUSA0 (with
   c = GMWC_MINUSA0 or c = GMWC_A3 the states obtained from different seeds
   are arithmetically related), and to seed a SplitMix generator and use
   its output to fill x, y, and z: do not use related seeds (e.g.,
   consecutive integers) directly, as they generate correlated sequences.
   To obtain different streams, you can also use the jump functions. */
uint64_t x, y, z, c;

uint64_t next(void) {
	const __uint128_t t = GMWC_A3 * (__uint128_t)x + c;
	x = y;
	y = z;
	z = GMWC_A0INV * (uint64_t)t;
	c = (t + GMWC_MINUSA0 * (__uint128_t)z) >> 64;
	return z ^ (uint64_t)(t >> 64);
}


/* The following jump functions use a minimal multiprecision library. They
   map the state to the integer
   h = c * b^3 - GMWC_MINUSA0 * (x + y * b + z * b^2) mod m, where b = 2^64:
   since each call to next() multiplies h modulo m by A = b^-1 mod m,
   jumping ahead by n steps multiplies h by A^n mod m. The state is then
   recovered as x + y * b + z * b^2 = rec * h mod b^3, where
   rec = (-GMWC_MINUSA0)^-1 mod b^3, and
   c = (h + GMWC_MINUSA0 * (x + y * b + z * b^2)) / b^3. */

#define MP_SIZE 5
#include "mp.c"

static uint64_t mod[MP_SIZE] = { GMWC_MINUSA0, 0, 0, GMWC_A3 };
static uint64_t minusa0[MP_SIZE] = { GMWC_MINUSA0 };
// (-GMWC_MINUSA0)^-1 mod 2^192. This is used to recover x, y, z.
static uint64_t rec[MP_SIZE] = { 0xac8f71df86d6542d, 0x1d6eb5a8e2770044, 0xab086186e8e70614 };


/* This is the jump function for the generator. It is equivalent
   to 2^128 calls to next(); it can be used to generate 2^127
   non-overlapping subsequences for parallel computations.

   Equivalent C++ Boost multiprecision code:

   cpp_int b = cpp_int(1) << 64;
   cpp_int b3 = b * b * b;
   cpp_int m = GMWC_A3 * b3 + GMWC_MINUSA0;
   cpp_int r = cpp_int("0xf5b5310ac371f0da4f38f128595669a797a4e5a3d3c9847bdb4a76eba65a977c");
   cpp_int s = ((m + c * b3 - GMWC_MINUSA0 * (x + y * b + z * b * b)) * r) % m;

   cpp_int xyz = s * cpp_int("0xab086186e8e706141d6eb5a8e2770044ac8f71df86d6542d");
   x = uint64_t(xyz);
   y = uint64_t(xyz >> 64);
   z = uint64_t(xyz >> 128);
   c = uint64_t((s + (xyz % b3) * GMWC_MINUSA0) >> 192);
*/

void jump(void) {
	static uint64_t jump[MP_SIZE] = { 0xdb4a76eba65a977c, 0x97a4e5a3d3c9847b, 0x4f38f128595669a7, 0xf5b5310ac371f0da };

	uint64_t state[MP_SIZE] = { 0, 0, 0, c };
	uint64_t xyz[MP_SIZE] = { x, y, z };
	mp_mul(xyz, minusa0, NULL);
	mp_sub(state, xyz, mod);
	mp_mul(state, jump, mod);
	memcpy(xyz, state, sizeof xyz);
	mp_mul(xyz, rec, NULL);

	x = xyz[0];
	y = xyz[1];
	z = xyz[2];
	xyz[3] = xyz[4] = 0;
	mp_mul(xyz, minusa0, NULL);
	mp_add(state, xyz, NULL);
	c = state[3];
}


/* This is the long-jump function for the generator. It is equivalent to
   2^192 calls to next(); it can be used to generate 2^63 starting points,
   from each of which jump() will generate 2^64 non-overlapping
   subsequences for parallel distributed computations.

   Equivalent C++ Boost multiprecision code:

   cpp_int b = cpp_int(1) << 64;
   cpp_int b3 = b * b * b;
   cpp_int m = GMWC_A3 * b3 + GMWC_MINUSA0;
   cpp_int r = cpp_int("0x9451fc15e8072ae4cc2b1b529d0106541456221459220092be888462fd857d89");
   cpp_int s = ((m + c * b3 - GMWC_MINUSA0 * (x + y * b + z * b * b)) * r) % m;

   cpp_int xyz = s * cpp_int("0xab086186e8e706141d6eb5a8e2770044ac8f71df86d6542d");
   x = uint64_t(xyz);
   y = uint64_t(xyz >> 64);
   z = uint64_t(xyz >> 128);
   c = uint64_t((s + (xyz % b3) * GMWC_MINUSA0) >> 192);
*/

void long_jump(void) {
	static uint64_t long_jump[MP_SIZE] = { 0xbe888462fd857d89, 0x1456221459220092, 0xcc2b1b529d010654, 0x9451fc15e8072ae4 };

	uint64_t state[MP_SIZE] = { 0, 0, 0, c };
	uint64_t xyz[MP_SIZE] = { x, y, z };
	mp_mul(xyz, minusa0, NULL);
	mp_sub(state, xyz, mod);
	mp_mul(state, long_jump, mod);
	memcpy(xyz, state, sizeof xyz);
	mp_mul(xyz, rec, NULL);

	x = xyz[0];
	y = xyz[1];
	z = xyz[2];
	xyz[3] = xyz[4] = 0;
	mp_mul(xyz, minusa0, NULL);
	mp_add(state, xyz, NULL);
	c = state[3];
}


/* The following arbitrary-jump functions compute A^n mod m by modular
   exponentiation. Since m has no special form, A is given explicitly. */

static uint64_t base[MP_SIZE] = { 0x0023a093509657d4, 0, 0x38472a068a14644e, 0xac330795a8035308 };   // A = b^-1 mod m

static void jumpmul(const uint64_t * const r) {
	uint64_t state[MP_SIZE] = { 0, 0, 0, c };
	uint64_t xyz[MP_SIZE] = { x, y, z };
	mp_mul(xyz, minusa0, NULL);
	mp_sub(state, xyz, mod);
	mp_mul(state, r, mod);
	memcpy(xyz, state, sizeof xyz);
	mp_mul(xyz, rec, NULL);

	x = xyz[0];
	y = xyz[1];
	z = xyz[2];
	xyz[3] = xyz[4] = 0;
	mp_mul(xyz, minusa0, NULL);
	mp_add(state, xyz, NULL);
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
   cpp_int m = GMWC_A3 * b3 + GMWC_MINUSA0;
   cpp_int A("0xac330795a803530838472a068a14644e00000000000000000023a093509657d4");  // = b^-1 mod m
   cpp_int r = powm(A, cpp_int(c) << e, m);                                          // A^(c * 2^e) mod m
   cpp_int s = ((m + carry * b3 - GMWC_MINUSA0 * (x + y * b + z * b * b)) * r) % m;

   cpp_int xyz = s * cpp_int("0xab086186e8e706141d6eb5a8e2770044ac8f71df86d6542d");
   x = uint64_t(xyz);
   y = uint64_t(xyz >> 64);
   z = uint64_t(xyz >> 128);
   carry = uint64_t((s + (xyz % b3) * GMWC_MINUSA0) >> 192);
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
   cpp_int m = GMWC_A3 * b3 + GMWC_MINUSA0;
   cpp_int A("0xac330795a803530838472a068a14644e00000000000000000023a093509657d4");  // = b^-1 mod m
   cpp_int r = powm(A, n, m);                                                        // A^n mod m
   cpp_int s = ((m + carry * b3 - GMWC_MINUSA0 * (x + y * b + z * b * b)) * r) % m;

   cpp_int xyz = s * cpp_int("0xab086186e8e706141d6eb5a8e2770044ac8f71df86d6542d");
   x = uint64_t(xyz);
   y = uint64_t(xyz >> 64);
   z = uint64_t(xyz >> 128);
   carry = uint64_t((s + (xyz % b3) * GMWC_MINUSA0) >> 192);
*/

void jump_n(const uint64_t jump[MP_SIZE - 1]) {
	uint64_t r[MP_SIZE];
	mp_powmod_arr(r, base, jump, MP_SIZE - 1, mod);
	jumpmul(r);
}
