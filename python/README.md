Code for the Selection of Multipliers
=====================================

The code in this directory selects multipliers from the databases of
candidates. It uses [NumPy](https://numpy.org/) and
[pandas](https://pandas.pydata.org/).

- `selMWC.py N <DATABASE` prints the `N` candidates with the largest weighted
  geometric mean *G*₈ of the figures of merit *f*₂, …, *f*₈, with weights
  1/(*t* − 1) (those of the harmonic score of Steele and Vigna), followed by
  *G*₈. The ranking does not depend on the normalization of the figures of
  merit, nor on the other candidates. For GMWC databases, use the symbolic
  link `selGMWC.py`: the script chooses the layout of the database from its
  name.
