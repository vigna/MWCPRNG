/* Test of mul_add.h in CUDA device code: compares mul_add() computed on the
   GPU with __uint128_t arithmetic on the host on all triples of special
   values, on random triples, and on the steps of the MWC generators with
   64-bit digits. Compile from the parent directory with

   nvcc -O3 -arch=sm_XY -I. test/mul_add.cu

   (-arch=native from CUDA 11.6 on), adding -DMUL_ADD_NO_INT128 or
   -DMUL_ADD_PORTABLE to test the other implementations. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "mul_add.h"

#define CHECK(call) do { const cudaError_t e = (call); if (e != cudaSuccess) { printf("%s: %s\n", #call, cudaGetErrorString(e)); return 1; } } while (0)

static const char * const names[] = { "unsigned __int128", "__umul64hi()", "schoolbook" };

__device__ int method(void) {
#if defined(MUL_ADD_INT128)
	return 0;
#elif defined(MUL_ADD_UMUL64HI)
	return 1;
#else
	return 2;
#endif
}

__global__ void triples(const uint64_t *a, const uint64_t *x, const uint64_t *c, uint64_t *lo, uint64_t *hi, const int n, int *m) {
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i == 0) *m = method();
	if (i < n) lo[i] = mul_add(a[i], x[i], c[i], &hi[i]);
}

/* Each thread runs n steps of the MWC generator with lag k and multiplier a
   on its state, made of k digits followed by the carry. */
__global__ void steps(const uint64_t a, const int k, const long n, uint64_t *state, const int threads) {
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= threads) return;
	uint64_t * const s = state + i * (k + 1), d[3];
	for (int j = 0; j < k; j++) d[j] = s[j];
	uint64_t c = s[k];
	for (long t = 0; t < n; t++) {
		uint64_t hi;
		const uint64_t lo = mul_add(a, d[0], c, &hi);
		for (int j = 0; j < k - 1; j++) d[j] = d[j + 1];
		d[k - 1] = lo;
		c = hi;
	}
	for (int j = 0; j < k; j++) s[j] = d[j];
	s[k] = c;
}

static uint64_t z = 0x9e3779b97f4a7c15;

static uint64_t splitmix64(void) {
	uint64_t r = (z += 0x9e3779b97f4a7c15);
	r = (r ^ (r >> 30)) * 0xbf58476d1ce4e5b9;
	r = (r ^ (r >> 27)) * 0x94d049bb133111eb;
	return r ^ (r >> 31);
}

int main(void) {
	static const uint64_t v[] = { 0, 1, 2, 0xffffffff, 0x100000000, 0x100000001, 0x7fffffffffffffff, 0x8000000000000000,
		0xfffffffeffffffff, 0xffffffff00000000, 0xfffffffffffffffe, 0xffffffffffffffff,
		0xff6da236274e252e, 0xff4e24dd418f5433, 0xff32e7e2bdcb0a33 };
	const int nv = sizeof v / sizeof *v, nr = 1 << 22, n = nv * nv * nv + nr;
	uint64_t *a, *x, *c, *lo, *hi;
	int *m;
	CHECK(cudaMallocManaged(&a, n * sizeof *a));
	CHECK(cudaMallocManaged(&x, n * sizeof *x));
	CHECK(cudaMallocManaged(&c, n * sizeof *c));
	CHECK(cudaMallocManaged(&lo, n * sizeof *lo));
	CHECK(cudaMallocManaged(&hi, n * sizeof *hi));
	CHECK(cudaMallocManaged(&m, sizeof *m));
	int p = 0;
	for (int i = 0; i < nv; i++)
		for (int j = 0; j < nv; j++)
			for (int k = 0; k < nv; k++) a[p] = v[i], x[p] = v[j], c[p++] = v[k];
	for (int i = 0; i < nr; i++) {
		a[p] = splitmix64();
		x[p] = i % 3 == 1 ? v[i % nv] : splitmix64();
		c[p++] = i % 3 == 2 ? v[i % nv] : splitmix64();
	}
	triples<<<(n + 255) / 256, 256>>>(a, x, c, lo, hi, n, m);
	CHECK(cudaGetLastError());
	CHECK(cudaDeviceSynchronize());
	long errors = 0;
	for (int i = 0; i < n; i++) {
		const __uint128_t t = (__uint128_t)a[i] * x[i] + c[i];
		if (lo[i] != (uint64_t)t || hi[i] != (uint64_t)(t >> 64)) {
			if (errors++ < 10) printf("mul_add(%016" PRIx64 ", %016" PRIx64 ", %016" PRIx64 ") = %016" PRIx64 "%016" PRIx64 "\n", a[i], x[i], c[i], hi[i], lo[i]);
		}
	}

	static const uint64_t mult[] = { 0xff6da236274e252e, 0xff4e24dd418f5433, 0xff32e7e2bdcb0a33 };
	const int threads = 1024;
	const long ns = 100000;
	uint64_t *s, *e = (uint64_t *)malloc(threads * 4 * sizeof *e);
	CHECK(cudaMallocManaged(&s, threads * 4 * sizeof *s));
	for (int k = 1; k <= 3; k++) {
		const uint64_t A = mult[k - 1];
		for (int i = 0; i < threads; i++) {
			for (int j = 0; j < k; j++) s[i * (k + 1) + j] = splitmix64();
			s[i * (k + 1) + k] = 1 + ((__uint128_t)splitmix64() * (A - 2) >> 64);
		}
		memcpy(e, s, threads * (k + 1) * sizeof *e);
		steps<<<threads / 256, 256>>>(A, k, ns, s, threads);
		CHECK(cudaGetLastError());
		CHECK(cudaDeviceSynchronize());
		for (int i = 0; i < threads; i++) {
			uint64_t * const t = e + i * (k + 1);
			for (long j = 0; j < ns; j++) {
				const __uint128_t r = (__uint128_t)A * t[0] + t[k];
				for (int l = 0; l < k - 1; l++) t[l] = t[l + 1];
				t[k - 1] = (uint64_t)r;
				t[k] = (uint64_t)(r >> 64);
			}
			for (int j = 0; j <= k; j++)
				if (s[i * (k + 1) + j] != t[j]) {
					if (errors++ < 10) printf("thread %d of the generator with multiplier %016" PRIx64 " differs\n", i, A);
					break;
				}
		}
	}
	printf("mul_add on the device (%s): %s\n", names[*m], errors ? "FAIL" : "OK");
	return errors != 0;
}
