#!/usr/bin/python3

#  Written in 2021 by Sebastiano Vigna
#
# To the extent possible under law, the author has dedicated all copyright
# and related and neighboring rights to this software to the public domain
# worldwide. This software is distributed without any warranty.
#
# See <http://creativecommons.org/publicdomain/zero/1.0/>.

# Prints the ITEMS candidates with the largest weighted geometric mean of the
# figures of merit f₂, …, f₈, with weights 1/(t − 1) (those of the harmonic score
# of Steele and Vigna). The ranking does not depend on the normalization of the
# figures of merit, nor on the other candidates.

import sys
import os
import numpy as np
import pandas as pd

if len(sys.argv) < 2:
    sys.stderr.write("Usage: %s ITEMS\n" % sys.argv[0])
    exit()

name = os.path.basename(__file__)

if name == "selMWC.py":
    mults = ['a']
elif name == "selGMWC.py":
    mults = ['-a₀', 'a₀⁻¹', 'a']
else:
    raise Exception("This script must be called selMWC.py or selGMWC.py")

labels = ['f₂', 'f₃', 'f₄', 'f₅', 'f₆', 'f₇', 'f₈']
data = pd.read_csv(sys.stdin, sep='\t', names=['M'] + mults + labels, index_col=False,
                   usecols=range(1 + len(mults) + len(labels)), dtype={c: 'str' for c in ['M'] + mults},
                   on_bad_lines='skip')
for c in labels:
    data[c] = pd.to_numeric(data[c], errors='coerce')
data = data[data['M'].str.startswith('0x', na=False) & (data[labels] > 0).all(axis=1)]

# Candidates are drawn with replacement, so some might appear twice
data = data.drop_duplicates(subset='M')
weights = 1 / np.arange(1, len(labels) + 1)   # 1/(t - 1), t = 2..8
data = data.assign(G=np.exp(np.log(data[labels]) @ weights / weights.sum())).nlargest(int(sys.argv[1]), 'G')

for index, row in data.iterrows():
    print("\t".join([row['M']] + [row[c] for c in mults] + ["%g" % row[c] for c in labels + ['G']]))
