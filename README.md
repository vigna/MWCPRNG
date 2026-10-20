# MWCPRNG

This repository contains data and code associated with the paper
“Scrambled Multiply-with-Carry Pseudorandom Number Generators”, by Sebastiano
Vigna.

The directory `search` contains the code used to sample candidate
multipliers and to compute their spectral figures of merit.

The directory `python` contains the Python code used to select the
multipliers reported in the paper.

The directory `prngs` contains the C code of the generators suggested in
the paper, which is also available at <https://prng.di.unimi.it/>, and a
test suite (`prngs/test/check.sh`) comparing their outputs and jumps with an
independent Python model.

The databases of candidate multipliers are available at
<https://vigna.di.unimi.it/MWC.tar.zst>. (The size of the zstd file is
approximately 1.8 gigabytes; the uncompressed data is approximately 7
gigabytes.)

The databases contain, for each type, state size, and base, the candidates
on which the spectral test has been run (e.g., `MWC-128[64]` for MWC
generators with 128 bits of state, carry included, and base 2⁶⁴, or
`GMWC-64[32]` for GMWC generators with 64 bits of state and base 2³²). For
each candidate we list, tab-separated, the modulus *m* of the associated
multiplicative congruential generator, the multiplier *a* (for GMWC
generators, *e* = −*a*₀, *a*₀⁻¹ mod *b*, and *a* = *a*₁), all in hexadecimal
form, and finally the figures of merit from dimension 2 to dimension 8.

For `MWC-64[32]` and `MWC-128[32]` the search is exhaustive over the
multipliers in [2³² − 2²⁴ . . 2³²); in all other cases candidates are random,
so some might appear twice.

The figures of merit have been computed using
[LatticeTester](https://github.com/umontreal-simul/latticetester), which
prints six significant digits.
