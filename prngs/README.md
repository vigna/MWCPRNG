Code of the Suggested Generators
================================

The code in this directory implements the generators suggested in the
paper, which are also available at <https://prng.di.unimi.it/>. It needs a C
compiler supporting 128-bit integers (`__uint128_t`).

All generators are scrambled by X: each step computes
*t* = *a* *x*ₙ₋ₖ + *c*ₙ₋₁ (for GMWC generators, *a* = *a*₁), from which the new
digit *x*ₙ and the new carry *c*ₙ are derived, and X returns
*x*ₙ ⊕ ⌊*t* / *b*⌋, the new digit xored with the upper half of *t*. For MWC
generators ⌊*t* / *b*⌋ = *c*ₙ, so X reduces to the scrambler of David
Thomas's MWC64X generator, and the output can be computed on the current
state.

| Generator | *b* | *k* | Constants | Jump | Long jump |
|---|---|---|---|---|---|
| `MWC128X.c` | 2⁶⁴ | 1 | `MWC_A1 = 0xff6da236274e252e` | 2⁶⁴ | 2⁹⁶ |
| `MWC192X.c` | 2⁶⁴ | 2 | `MWC_A2 = 0xff4e24dd418f5433` | 2⁹⁶ | 2¹⁴⁴ |
| `MWC256X.c` | 2⁶⁴ | 3 | `MWC_A3 = 0xff32e7e2bdcb0a33` | 2¹²⁸ | 2¹⁹² |
| `GMWC128X.c` | 2⁶⁴ | 1 | `GMWC_MINUSA0 = 0x68935f57`, `GMWC_A1 = 0xfff0cfb4ed6366fa` | 2⁶⁴ | 2⁹⁶ |
| `GMWC256X.c` | 2⁶⁴ | 3 | `GMWC_MINUSA0 = 0x34dabaf938e45b`, `GMWC_A3 = 0xff76e600f25ae046` | 2¹²⁸ | 2¹⁹² |
| `MWC64X_32.c` | 2³² | 1 | `MWC_A1 = 0xfff145f8` | 2³² | 2⁴⁸ |
| `MWC128X_32.c` | 2³² | 3 | `MWC_A3 = 0xffcccc66` | 2⁶⁴ | 2⁹⁶ |

Besides `jump()` and `long_jump()`, every generator provides `jump_ce(c, e)`,
which is equivalent to *c* · 2ᵉ calls to `next()`, and `jump_n()`, which
jumps ahead by an arbitrary number of steps. The jump functions use the
minimal multiprecision library `mp.c`, except for those of `MWC64X_32.c`,
whose modulus fits in 64 bits, which use 128-bit arithmetic instead.

To seed a generator, set the carry to one (for GMWC generators, to
2*e*, where *e* = −*a*₀ = `GMWC_MINUSA0`) and fill the other words with the
output of a SplitMix generator.

`test/check.sh` compiles `test/test.c` with each generator at `-O0` and
`-O3`, checks that `jump_ce(1, 16)` is equivalent to 2¹⁶ calls to `next()` and
that `jump()` and `long_jump()` are equivalent to the corresponding calls to
`jump_ce()`, and compares the outputs and the results of all jump functions
with an independent Python model (`test/check.py`), which also checks that
the modulus *m* and (*m* − 1) / 2 are prime.
