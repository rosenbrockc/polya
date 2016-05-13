# Reproducing the Figures

Below we list the figures in the
[paper](https://github.com/rosenbrockc/polya/blob/master/docs/polyaenum.pdf)
that deal with timing experiments and include details on how to
generate the data for each one.

## Figure 5: Algorithm Scaling with Colors

We used a transitive group constructed from $S_3\times S_4$ (thus with
144 group elements). This group is stored in the file
[`fortran/tests/group.out.cr1`](https://github.com/rosenbrockc/polya/blob/master/fortran/tests/group.out.cr1). Additionally,
the generators are contained in the same folder in a file
`generators.in.cr1`.

| Colors |         Stochiometry         |
|:------:|:----------------------------:|
|    2   | 10 10                        |
|    3   | 10  5  5                     |
|    4   | 5  5  5  5                   |
|    5   | 4  4  4  4  4                |
|    6   | 4  4  4  4  2  2             |
|    7   | 4  4  4  2  2  2  2          |
|    8   | 4  4  2  2  2  2  2  2       |
|    9   | 4  2  2  2  2  2  2  2  2    |
|   10   | 2  2  2  2  2  2  2  2  2  2 |

The Polya value for the first entry in the table can be generated
using `./polya.out 10 10 -group group.out.cr1`.

## Figure 6: Algorithm Scaling with Number of Elements in the Set

For a fixed color selection (only 2 colors), we adjust the
stoichiometry as follows:

|  Group | Stochiometry |
|:------:|:------------:|
|   fg1  |    6     6   |
|   fg2  |    9     9   |
|   fg3  |   10    10   |
|   fg4  |   12    12   |
|   fg5  |   15    15   |
|   fg6  |   20    20   |
|   fg7  |   55    17   |
|   fg8  |  110    10   |     

The groups are located at
[`fortran/tests/`](https://github.com/rosenbrockc/polya/blob/master/fortran/tests/),
with files labeled `group.out.{group}`, with `{group}` being one of
the `fg*` in the table above. The groups are all isomorphic to
$S_3\times S_4$ (thus with 144 group elements), but constructed
carefully to be transitive on the sites. Additionally, the generators
are contained in the same folder in a file `generators.in.fg*`.

The Polya value for the first entry in the table can be generated
using `./polya.out 6 6 -group group.out.fg1`.

## Figure 7: Algorithm Scaling with Group Size

We used the unique permutation groups arising from all derivative
super structures of a simple cubic lattice for a given number of sites
in the unit cell. The groups can be found in the `group.out.sc*` files
in
[`fortran/tests/`](https://github.com/rosenbrockc/polya/blob/master/fortran/tests/). Additionally,
the generators are contained in the same folder in a file
`generators.in.sc_*`.

|  Group | Stochiometry |
|:------:|:------------:|
| sc_16* |    8     8   |
| sc_24* |   12    12   |
| sc_32* |   16    16   |
| sc_48* |   24    24   |
| sc_60* |   30    30   |