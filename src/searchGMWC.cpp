/*  Written in 2020-2021 by Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

/* Computes valid constants for a Goreski-Klapper GMWC. A few possible state
   and output bit size are available (see the list of #ifdef below). For each
   set of constants, the output is given by four positive numbers: m (the prime
   modulus of the induced LCG), -a_0, the modular inverse of a_0, and a_1. Note
   that m = a1 * b^k - a_0, where b is the base (i.e., the word size) and k the
   number of words of state (not including the carry). m and (m-1)/2 will be
   both primes.

   Beside a seed for the search, you can specify a number of bits that will be
   used to limit the size of a_0. This can be useful if a smaller a_0 results
   in a faster multiplication. However, the main advantage of GMWCs is that the
   binary representation of the prime m does not end with ...00001, so reducing
   the size of excessively a_0 will reduce the quality of the generator. */

#include <iostream>
#include <NTL/LLL.h>

using namespace NTL;
using namespace std;

static inline uint64_t rotl(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}

static uint64_t s[4] = { -1ULL, 1, -1ULL, 1 };

static uint64_t next(void) {
	const uint64_t result = rotl(s[0] + s[3], 23) + s[0];

	const uint64_t t = s[1] << 17;

	s[2] ^= s[0];
	s[3] ^= s[1];
	s[1] ^= s[2];
	s[0] ^= s[3];

	s[2] ^= t;

	s[3] = rotl(s[3], 45);

	return result;
}

int main(int argc, char *argv[]) {
#if defined(BITS_32_64)
	if (argc < 2) {
		fprintf(stderr, "Usage: %s SEED [BITS]\n", argv[0]);
		return 1;
	}

	s[0] = strtoull(argv[1], NULL, 0);

	const int bits = argc == 2 ? 32 : atoi(argv[2]);
	if (bits <= 0 || bits > 32)  {
		fprintf(stderr, "Invalid number of bits: %d\n", bits);
		return 1;
	}
	const uint64_t mask = UINT64_C(-1) >> 64 - bits;

	ZZ b = conv<ZZ>(1) << 32;
	ZZ m, minus_a0, a1, a0_inv;

	for(;;) {
		a1 = conv<ZZ>((next() & 0xFFFFFFFF) | 0xFF000000);
		minus_a0 = conv<ZZ>((next() & mask) | 1);
		if (a1 + minus_a0 >= b) continue;
		m = a1 * b + minus_a0;
		if (ProbPrime((m - 1) / 2, 2) && ProbPrime(m, 100) && ProbPrime((m - 1) / 2, 100)) {
			cout << m << "\t" << minus_a0 << "\t" << InvMod(b - minus_a0, b) << "\t" << a1 << endl;
		}
	}
#elif defined(BITS_32_128)
	if (argc < 2) {
		fprintf(stderr, "Usage: %s SEED [BITS]\n", argv[0]);
		return 1;
	}

	s[0] = strtoull(argv[1], NULL, 0);

	const int bits = argc == 2 ? 32 : atoi(argv[2]);
	if (bits <= 0 || bits > 32)  {
		fprintf(stderr, "Invalid number of bits: %d\n", bits);
		return 1;
	}
	const uint64_t mask = UINT64_C(-1) >> 64 - bits;

	ZZ b = conv<ZZ>(1) << 32;
	ZZ bbb = b * b * b;
	ZZ m, minus_a0, a1, a0_inv;

	for(;;) {
		a1 = conv<ZZ>((next() & 0xFFFFFFFF) | 0xFF000000);
		minus_a0 = conv<ZZ>((next() & mask) | 1);
		if (a1 + minus_a0 >= b) continue;
		m = a1 * bbb + minus_a0;
		if (ProbPrime((m - 1) / 2, 2) && ProbPrime(m, 100) && ProbPrime((m - 1) / 2, 100)) {
			cout << m << "\t" << minus_a0 << "\t" << InvMod(b - minus_a0, b) << "\t" << a1 << endl;
		}
	}
#elif defined(BITS_64_128)
	if (argc < 2) {
		fprintf(stderr, "Usage: %s SEED [BITS]\n", argv[0]);
		return 1;
	}

	s[0] = strtoull(argv[1], NULL, 0);

	const int bits = argc == 2 ? 64 : atoi(argv[2]);
	if (bits <= 0 || bits > 64)  {
		fprintf(stderr, "Invalid number of bits: %d\n", bits);
		return 1;
	}
	const uint64_t mask = UINT64_C(-1) >> 64 - bits;

	ZZ b = conv<ZZ>(1) << 64;
	ZZ m, minus_a0, a1, a0_inv;

	for(;;) {
		a1 = conv<ZZ>((next() & 0xFFFFFFFFFFFFFFFF) | 0xFF00000000000000);
		minus_a0 = conv<ZZ>((next() & mask) | 1);
		if (a1 + minus_a0 >= b) continue;
		m = a1 * b + minus_a0;
		if (ProbPrime((m - 1) / 2, 2) && ProbPrime(m, 100) && ProbPrime((m - 1) / 2, 100)) {
			cout << m << "\t" << minus_a0 << "\t" << InvMod(b - minus_a0, b) << "\t" << a1 << endl;
		}
	}
#elif defined(BITS_64_256)
	if (argc < 2) {
		fprintf(stderr, "Usage: %s SEED [BITS]\n", argv[0]);
		return 1;
	}

	s[0] = strtoull(argv[1], NULL, 0);

	const int bits = argc == 2 ? 64 : atoi(argv[2]);
	if (bits <= 0 || bits > 64)  {
		fprintf(stderr, "Invalid number of bits: %d\n", bits);
		return 1;
	}
	const uint64_t mask = UINT64_C(-1) >> 64 - bits;

	ZZ b = conv<ZZ>(1) << 64;
	ZZ bbb = b * b * b;
	ZZ m, minus_a0, a1, a0_inv;

	for(;;) {
		a1 = conv<ZZ>((next() & 0xFFFFFFFFFFFFFFFF) | 0xFF00000000000000);
		minus_a0 = conv<ZZ>((next() & mask) | 1);
		if (a1 + minus_a0 >= b) continue;
		m = a1 * bbb + minus_a0;
		if (ProbPrime((m - 1) / 2, 2) && ProbPrime(m, 100) && ProbPrime((m - 1) / 2, 100)) {
			cout << m << "\t" << minus_a0 << "\t" << InvMod(b - minus_a0, b) << "\t" << a1 << endl;
		}
	}
#elif defined(BITS_8_32)
	ZZ b = conv<ZZ>(1) << 8;
	ZZ bbb = b * b * b;
	ZZ m, minus_a0, a1, a0_inv;

	for(int a = 0xC0; a <= 0xFF; a++) {
		a1 = conv<ZZ>(a);
		for(int t = 1; t <= 0xFF - a; t++) {
			minus_a0 = conv<ZZ>(t);
			m = a1 * bbb + minus_a0;
			if (ProbPrime(m, 100) && ProbPrime((m - 1) / 2, 100)) {
				cout << m << "\t" << minus_a0 << "\t" << InvMod(b - minus_a0, b) << "\t" << a1 << endl;
			}
		}
	}
#elif defined(BITS_16_32)
	ZZ b = conv<ZZ>(1) << 16;
	ZZ m, minus_a0, a1, a0_inv;

	for(int a = 0xF000; a <= 0xFFFF; a++) {
		a1 = conv<ZZ>(a);
		for(int t = 1; t <= 0xFFFF - a; t++) {
			minus_a0 = conv<ZZ>(t);
			m = a1 * b + minus_a0;
			if (ProbPrime(m, 100) && ProbPrime((m - 1) / 2, 100)) {
				cout << m << "\t" << minus_a0 << "\t" << InvMod(b - minus_a0, b) << "\t" << a1 << endl;
			}
		}
	}
#else
#error "Please define the state and output number of bits"
#endif

	return 0;
}
