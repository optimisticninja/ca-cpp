# ca

*In development*

Playground for finding reversible cellular automata transforms in a 1D N-neighborhood cellular space.

Currently implemented:

* Elementary CA
* Game of Life

TODO:

* Game of life as base for 2D
* 1D/2D margolus, critters, tron rule for block implementation
    - See BlockCA1D in ca1d.h, partitioning and evolution loop written - need to implement rule
* [Dynamic Rules in Elementary Celluar Automata](https://github.com/gojakuch/dynamic-rule-cellular-automata/blob/main/paper_SOECA.pdf)
    - TODO: Create second-order rules from this implementation with similar "rule compiler"
* Devise rule seeding/using CA PRNG (src/ca/capring.cc) for translation

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
* [A Study and Comparison of First and Second Order Cellular Automata with Examples](https://ntnuopen.ntnu.no/ntnu-xmlui/bitstream/handle/11250/258721/351877_FULLTEXT01.pdf?sequence=3)
* [Reversible CA](http://www.jcasim.de/main/node10.html)
