# ca

*In development*

Playground for finding reversible cellular automata transforms in a 1D N-neighborhood cellular space.

Currently implemented:

* Elementary CA

## Terminology

There is a lot of varying terminology between papers, this is my map for consistent code and more generic terms in some cases

| term | definition | type |
| --- | --- | --- |
| gateway key | the encoding/configuration of the cellular automata (initial state, partition size, total cells, boundary rule, etc.) |
| partition | the operating area (including the cell) when the CA is sequentially (over a block or neighborhood) updating state |
| partition size | the total number of cells in the partition  |
| neighborhood | the cells surrounding target cell in the operating area | CA_1D
| block | a non-overlapping partition of the CA state to operate on. instead of updating at the cellular level, the whole block is updated | CA_1D_BLOCK
| bias | when dealing with even number partitions, this is the side to pull the extra cell from (right or left)                        |       CA_1D*                                                          |

## References

* [Cellular Automata Transforms: Theory and Applications in Multimedia Compression, Encryption, and Modeling](https://www.amazon.com/Cellular-Automata-Transforms-Applications-Compression/dp/0792378571)
* [Pseudorandom Pattern Generation by a 4-Neighborhood Cellular Automata Based on a Probabilistic Analysis](http://www.iaeng.org/publication/IMECS2008/IMECS2008_pp1908-1913.pdf)
