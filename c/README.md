# MWCX and GMWCX generators, version 1.0

`MWC128X.c`, `MWC192X.c`, `MWC256X.c`, `GMWC128X.c`, and `GMWC256X.c` (version
1.0) replace the generators distributed at <https://prng.di.unimi.it/>
(`MWC128.c`, `MWC192.c`, `MWC256.c`, `GMWC128.c`, `GMWC256.c`, version 1.0);
`mp.c` is unchanged, except for three typos in the comments. `MWC128X_32.c` and
`MWC64X_32.c` are new: they are the generators with 32-bit digits and outputs (b
= 2^32) suggested in the paper, with lag 3 and 1, multipliers `0xffcccc66` and
`0xfff145f8`, and the scrambler X, which for them is necessary (without it, they
fail three-dimensional birthday-spacings tests in a few seconds). `MWC64X_32.c`
is Thomas's MWC64X with a different multiplier; since its modulus fits in 64
bits, its jump functions use 128-bit arithmetic instead of `mp.c`. With respect
to the published versions (downloaded on 16 October 2026):

- The multipliers (and, for GMWC, the parameter `GMWC_MINUSA0` and its inverse
  `GMWC_A0INV`) are those suggested in the paper (best weighted geometric mean
  of the figures of merit in dimensions 2 to 8):

  | Generator | Constants |
  |---|---|
  | `MWC128X` | `MWC_A1 = 0xff6da236274e252e` |
  | `MWC192X` | `MWC_A2 = 0xff4e24dd418f5433` |
  | `MWC256X` | `MWC_A3 = 0xff32e7e2bdcb0a33` |
  | `GMWC128X` | `GMWC_MINUSA0 = 0x68935f57`, `GMWC_A1 = 0xfff0cfb4ed6366fa` |
  | `GMWC256X` | `GMWC_MINUSA0 = 0x34dabaf938e45b`, `GMWC_A3 = 0xff76e600f25ae046` |

- The output is scrambled by X, a generalization of the scrambler of Thomas's
  MWC64X: the next digit xored with the upper half of the value `t` = a * x + c
  from which the next digit and carry are computed. For MWC generators, the
  upper half of `t` is the next carry, so X degenerates to Thomas's scrambler,
  and the output is computed, equivalently, as the most recent digit xored with
  the carry, both of the current state (as before, MWC generators compute the
  output on the current state). For GMWC generators, the digit xored with the
  carry would be biased, as the range of the carry depends on the digit (for
  `GMWC256X`, the most significant bit would be one with probability about
  1/2 - 2 * 10^-4).
- All constants depending on the parameters have been recomputed: the modulus in
  the comments, the jump and long-jump multipliers, the base `A = b^-1 mod m` of
  the GMWC arbitrary jumps, and the inverse of `-GMWC_MINUSA0` modulo 2^192 of
  `GMWC256X`, both in the code and in the equivalent Boost code.
- Comments: version 1.0, description of the scrambler, of the dependencies
  between consecutive steps, of the lattice defect and of the tests on it (as
  passes and failures), and of the choice of the multiplier; explanation of the
  jump functions in terms of the state integer; seeding advice from the paper
  (carry 1 for MWC, 2 * `GMWC_MINUSA0` for GMWC, other words from SplitMix;
  related seeds, such as consecutive integers, should be passed through
  SplitMix); the numbers of subsequences and starting points of the jump
  functions have been corrected (e.g., the long jump of `MWC192X` gives 2^47
  starting points, not 2^96).
- `next()` is no longer declared `inline` (an `inline` definition without an
  external one does not link at `-O0`), and `<string.h>` replaces `<strings.h>`
  (`GMWC256X.c` uses `memcpy()`).

`test/check.sh` compiles `test/test.c` with each generator at `-O0` and `-O3`,
checks that `jump_ce(1, 16)` is equivalent to 2^16 calls to `next()` and that
`jump()` and `long_jump()` are equivalent to the corresponding calls to
`jump_ce()`, and compares the outputs and the results of all jump functions with
the independent Python model in `test/check.py`, which also checks that each
modulus `m` and `(m - 1) / 2` are prime.
