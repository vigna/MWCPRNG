/* Test driver: compile with -DGENFILE="\"MWC128X.c\"" -DK=1 (lag), -DGMWC for
   GMWC generators, and -DW=32 for generators with 32-bit digits. */
#include <stdio.h>
#include <inttypes.h>
#include GENFILE

#ifndef W
#define W 64
#endif
#if W == 32
typedef uint32_t digit;
#else
typedef uint64_t digit;
#endif

#define D(v) (W / 4), (uint64_t)(v)
#if K == 1
#define SAVE(s) do { s[0] = x; s[1] = c; } while (0)
#define PRINT(tag) printf("%s %0*" PRIx64 " %0*" PRIx64 "\n", tag, D(x), D(c))
#elif K == 2
#define SAVE(s) do { s[0] = x; s[1] = y; s[2] = c; } while (0)
#define PRINT(tag) printf("%s %0*" PRIx64 " %0*" PRIx64 " %0*" PRIx64 "\n", tag, D(x), D(y), D(c))
#else
#define SAVE(s) do { s[0] = x; s[1] = y; s[2] = z; s[3] = c; } while (0)
#define PRINT(tag) printf("%s %0*" PRIx64 " %0*" PRIx64 " %0*" PRIx64 " %0*" PRIx64 "\n", tag, D(x), D(y), D(z), D(c))
#endif

static void seed(void) {
	x = (digit)0x0123456789abcdef;
#if K >= 2
	y = (digit)0xfedcba9876543210;
#endif
#if K >= 3
	z = (digit)0x0f1e2d3c4b5a6978;
#endif
#ifdef GMWC
	c = 2 * (uint64_t)GMWC_MINUSA0;
#else
	c = 1;
#endif
}

int main(void) {
	uint64_t s[4], t[4];
	seed();
	for (int i = 0; i < 8; i++) printf("out %0*" PRIx64 "\n", D(next()));
	seed(); jump(); PRINT("jump");
	seed(); long_jump(); PRINT("long_jump");
	seed(); jump_ce(3, 70); PRINT("jump_ce_3_70");
	seed(); jump_ce(12345, 0); PRINT("jump_ce_12345_0");
#ifdef MP_SIZE
	seed(); { uint64_t n[MP_SIZE - 1] = { 0x9e3779b97f4a7c15, 0x0123456789abcdef }; jump_n(n); } PRINT("jump_n");
#else
	seed(); jump_n(0x9e3779b97f4a7c15); PRINT("jump_n");
#endif
	/* Internal consistency: jump_ce(1, 16) == 2^16 calls to next(). */
	seed(); jump_ce(1, 16); SAVE(s);
	seed(); for (int i = 0; i < 65536; i++) next(); SAVE(t);
	int ok = 1; for (int i = 0; i < K + 1; i++) ok &= s[i] == t[i];
	printf("jump_ce(1,16)==2^16 next(): %s\n", ok ? "OK" : "FAIL");
	/* jump() == jump_ce(1, JE) and long_jump() == jump_ce(1, LE) */
	seed(); jump(); SAVE(s); seed(); jump_ce(1, JE); SAVE(t);
	ok = 1; for (int i = 0; i < K + 1; i++) ok &= s[i] == t[i];
	printf("jump()==jump_ce(1,%d): %s\n", JE, ok ? "OK" : "FAIL");
	seed(); long_jump(); SAVE(s); seed(); jump_ce(1, LE); SAVE(t);
	ok = 1; for (int i = 0; i < K + 1; i++) ok &= s[i] == t[i];
	printf("long_jump()==jump_ce(1,%d): %s\n", LE, ok ? "OK" : "FAIL");
	return 0;
}
