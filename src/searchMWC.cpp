/*  Written in 2020-2021 by Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

/* Computes valid constants for an MWC. A few possible state and output bit
   size are available (see the list of #ifdef below). The output is given by
   two positive numbers: m (the prime modulus of the induced LCG) and a. Note
   that m = a * b^k - 1, where b is the base (i.e., the word size) and k the
   number of words of state (not including the carry). m and (m-1)/2 will be
   both primes. */

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
	if (argc != 3) {
		fprintf(stderr, "Usage: %s FROM TO\n", argv[0]);
		return 1;
	}

	const uint64_t from  = strtoull(argv[1], NULL, 0);
	const uint64_t to  = strtoull(argv[2], NULL, 0);

	ZZ b = (conv<ZZ>(1) << 32);
	ZZ m;

	for(uint64_t a = from; a < to; a++) {
		m = conv<ZZ>(a) * b - 1;
		if (ProbPrime(m, 2) && ProbPrime((m - 1) / 2, 2) && ProbPrime(m, 100) && ProbPrime((m - 1) / 2, 100)) {
			cout << m << "\t" << a << endl;
		}
	}
#elif defined(BITS_32_128)
	if (argc != 3) {
		fprintf(stderr, "Usage: %s FROM TO\n", argv[0]);
		return 1;
	}

	const uint64_t from  = strtoull(argv[1], NULL, 0);
	const uint64_t to  = strtoull(argv[2], NULL, 0);

	ZZ b = (conv<ZZ>(1) << 32);
	ZZ bbb = b * b * b;
	ZZ m, a;

	for(uint64_t a = from; a < to; a++) {
		m = conv<ZZ>(a) * bbb - 1;
		if (ProbPrime((m - 1) / 2, 2) && ProbPrime(m, 100) && ProbPrime((m - 1) / 2, 100)) {
			cout << m << "\t" << a << endl;
		}
	}
#elif defined(BITS_64_128)
	if (argc != 2) {
		fprintf(stderr, "Usage: %s SEED\n", argv[0]);
		return 1;
	}

	s[0] = strtoull(argv[1], NULL, 0);
	ZZ b = (conv<ZZ>(1) << 64);
	ZZ m, a;

	for(;;) {
		a = conv<ZZ>((next() & 0xFFFFFFFFFFFFFFFF) | 0xFF00000000000000);
		m = a * b - 1;
		if (ProbPrime((m - 1) / 2, 2) && ProbPrime(m, 100) && ProbPrime((m - 1) / 2, 100)) {
			cout << m << "\t" << a << endl;
		}
	}
#elif defined(BITS_64_192)
	if (argc != 2) {
		fprintf(stderr, "Usage: %s SEED\n", argv[0]);
		return 1;
	}

	s[0] = strtoull(argv[1], NULL, 0);
	ZZ b = (conv<ZZ>(1) << 64);
	ZZ bb = b * b;
	ZZ m, a;

	for(;;) {
		a = conv<ZZ>((next() & 0xFFFFFFFFFFFFFFFF) | 0xFF00000000000000);
		m = a * bb - 1;
		if (ProbPrime((m - 1) / 2, 2) && ProbPrime(m, 100) && ProbPrime((m - 1) / 2, 100)) {
			cout << m << "\t" << a << endl;
		}
	}
#elif defined(BITS_64_256)
	if (argc != 2) {
		fprintf(stderr, "Usage: %s SEED\n", argv[0]);
		return 1;
	}

	s[0] = strtoull(argv[1], NULL, 0);
	ZZ b = (conv<ZZ>(1) << 64);
	ZZ bbb = b * b * b;
	ZZ m, a;

	for(;;) {
		a = conv<ZZ>((next() & 0xFFFFFFFFFFFFFFFF) | 0xFF00000000000000);
		m = a * bbb - 1;
		if (ProbPrime((m - 1) / 2, 2) && ProbPrime(m, 100) && ProbPrime((m - 1) / 2, 100)) {
			cout << m << "\t" << a << endl;
		}
	}
#elif defined(BITS_8_32)
	ZZ b = (conv<ZZ>(1) << 8);
	ZZ bbb = b * b * b;
	ZZ m;

	for(int a = 0xC0; a <= 0xFF; a++) {
		m = conv<ZZ>(a) * bbb - 1;
		if (ProbPrime(m, 100) && ProbPrime((m - 1) / 2, 100)) {
			cout << m << "\t" << a << endl;
		}
	}
#elif defined(BITS_16_32)
	ZZ b = (conv<ZZ>(1) << 16);
	ZZ m, a;

	for(int a = 0xF000; a <= 0xFFFF; a++) {
		m = conv<ZZ>(a) * b - 1;

		if (ProbPrime(m, 100) && ProbPrime((m - 1) / 2, 100)) {
			cout << m << "\t" << a << endl;
		}
	}
#elif defined(BITS_29_58)
	ZZ b = (conv<ZZ>(1) << 29);
	ZZ m, a;

	for(int a = 0x1FF00000; a <= 0x1FFFFFFF; a++) {
		m = a * b - 1;
		if (ProbPrime((m - 1) / 2, 2) && ProbPrime(m, 100) && ProbPrime((m - 1) / 2, 100)) {
			cout << m << "\t" << a << endl;
		}
	}
#elif defined(BITS_32_59)
	ZZ b = (conv<ZZ>(1) << 32);
	ZZ m, a;

	for(int a = 0x7F00000; a <= 0x7FFFFFF; a++) {
		m = a * b - 1;
		if (ProbPrime((m - 1) / 2, 2) && ProbPrime(m, 100) && ProbPrime((m - 1) / 2, 100)) {
			cout << m << "\t" << a << endl;
		}
	}
#elif defined(BITS_30_59)
	ZZ b = (conv<ZZ>(1) << 30);
	ZZ m, a;

	for(int a = 0x1FF00000; a <= 0x1FFFFFFF; a++) {
		m = a * b - 1;
		if (ProbPrime((m - 1) / 2, 2) && ProbPrime(m, 100) && ProbPrime((m - 1) / 2, 100)) {
			cout << m << "\t" << a << endl;
		}
	}
#elif defined(BITS_32_58)
	ZZ b = (conv<ZZ>(1) << 32);
	ZZ m, a;

	for(int a = 0x3F00000; a <= 0x3FFFFFF; a++) {
		m = a * b - 1;
		if (ProbPrime((m - 1) / 2, 2) && ProbPrime(m, 100) && ProbPrime((m - 1) / 2, 100)) {
			cout << m << "\t" << a << endl;
		}
	}
#else
#error "Please define the state and output number of bits"
#endif

	return 0;
}
