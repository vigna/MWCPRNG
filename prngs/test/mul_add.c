/* Test of mul_add.h: compares mul_add() with a reference implementation on
   all triples of special values, on random triples, and on the steps of the
   MWC generators with 64-bit digits. Compile with -DMUL_ADD_NO_INT128 or
   -DMUL_ADD_PORTABLE to test the other implementations. The reference uses
   __uint128_t if available; otherwise (e.g., with MSVC), a schoolbook
   multiplication written independently from that of mul_add.h. */
#include <stdio.h>
#include <inttypes.h>
#include "mul_add.h"

#if defined(MUL_ADD_INT128)
#define METHOD "128-bit integers"
#elif defined(MUL_ADD_UMUL128)
#define METHOD "_umul128()"
#elif defined(MUL_ADD_UMULH)
#define METHOD "__umulh()"
#else
#define METHOD "schoolbook"
#endif

static uint64_t s = 0x9e3779b97f4a7c15;

static uint64_t splitmix64(void) {
	uint64_t z = (s += 0x9e3779b97f4a7c15);
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
	return z ^ (z >> 31);
}

static void ref(const uint64_t a, const uint64_t x, const uint64_t c, uint64_t * const lo, uint64_t * const hi) {
#ifdef __SIZEOF_INT128__
	const __uint128_t t = (__uint128_t)a * x + c;
	*lo = (uint64_t)t;
	*hi = (uint64_t)(t >> 64);
#else
	const uint64_t a0 = (uint32_t)a, a1 = a >> 32, x0 = (uint32_t)x, x1 = x >> 32;
	const uint64_t p00 = a0 * x0, p01 = a0 * x1, p10 = a1 * x0, p11 = a1 * x1;
	/* (2^32 - 1) + (2^32 - 1) + (2^32 - 1)^2 = 2^64 - 1 */
	const uint64_t mid = (p00 >> 32) + (uint32_t)p10 + p01;
	const uint64_t l = (mid << 32 | (uint32_t)p00) + c;
	*lo = l;
	*hi = p11 + (p10 >> 32) + (mid >> 32) + (l < c);
#endif
}

static long errors;

static void check(const uint64_t a, const uint64_t x, const uint64_t c) {
	uint64_t hi, rlo, rhi;
	const uint64_t lo = mul_add(a, x, c, &hi);
	ref(a, x, c, &rlo, &rhi);
	if (lo != rlo || hi != rhi) {
		if (errors++ < 10) printf("mul_add(%016" PRIx64 ", %016" PRIx64 ", %016" PRIx64 ") = %016" PRIx64 "%016" PRIx64 "\n", a, x, c, hi, lo);
	}
}

/* Runs n steps of an MWC generator with lag k and multiplier a both with
   mul_add() and with the reference implementation. */
static void steps(const uint64_t a, const int k, const long n) {
	uint64_t d[3], e[3], c, f, r;
	ref(splitmix64(), a - 2, 0, &r, &c);
	f = ++c;
	for (int i = 0; i < k; i++) d[i] = e[i] = splitmix64();
	for (long i = 0; i < n; i++) {
		uint64_t hi, rlo, rhi;
		const uint64_t lo = mul_add(a, d[0], c, &hi);
		ref(a, e[0], f, &rlo, &rhi);
		for (int j = 0; j < k - 1; j++) d[j] = d[j + 1], e[j] = e[j + 1];
		d[k - 1] = lo;
		c = hi;
		e[k - 1] = rlo;
		f = rhi;
		if (d[k - 1] != e[k - 1] || c != f) {
			if (errors++ < 10) printf("step %ld of the generator with multiplier %016" PRIx64 " differs\n", i, a);
			return;
		}
	}
}

int main(void) {
	static const uint64_t v[] = { 0, 1, 2, 0xffffffff, 0x100000000, 0x100000001, 0x7fffffffffffffff, 0x8000000000000000,
		0xfffffffeffffffff, 0xffffffff00000000, 0xfffffffffffffffe, 0xffffffffffffffff,
		0xff6da236274e252e, 0xff4e24dd418f5433, 0xff32e7e2bdcb0a33 };
	const int n = (int)(sizeof v / sizeof *v);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++) check(v[i], v[j], v[k]);
	for (long i = 0; i < 10000000; i++) {
		const uint64_t a = splitmix64(), x = splitmix64(), c = splitmix64();
		check(a, x, c);
		check(a, x, v[i % n]);
		check(a, v[i % n], c);
	}
	steps(0xff6da236274e252e, 1, 10000000); /* MWC128X */
	steps(0xff4e24dd418f5433, 2, 10000000); /* MWC192X */
	steps(0xff32e7e2bdcb0a33, 3, 10000000); /* MWC256X */
	printf("mul_add (" METHOD "): %s\n", errors ? "FAIL" : "OK");
	return errors != 0;
}
