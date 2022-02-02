# ca-cpp

*In development*

Playground for different cellular automata

Currently implemented:

* 1D
	* Elementary CA
	* Second-order CA
* 2D
	* Game of Life

## Terminology

There is a lot of varying terminology between papers, this is my map for consistent code and more generic terms in some cases.

*NOTE*

Partitioning cellular automata (PCAs) will always be referred to as block cellular automata, as 'partition' is a good term for a low-level abstraction of the operating area during local transition (the block or neighborhood).

| term | definition | type |
| --- | --- | --- |
| gateway key | the encoding/configuration of the cellular automata (initial state, partition size, total cells, boundary rule, etc.) | * |
| partition | the operating area (sometimes including the cell) when the CA is sequentially (over a block or neighborhood) updating state | * |
| partition size | the total number of cells in the partition (block or neighborhood)  | * |
| neighborhood | the cells surrounding target cell in the operating area | CA_1D |
| block | a non-overlapping partition of the CA state to operate on. instead of updating at the cellular level, the whole block is updated | CA_1D_BLOCK |
| bias | when dealing with even number partitions, this is the side to pull the extra cell from (right or left) | CA_1D* |

## References

* [Cellular Automata Transforms: Theory and Applications in Multimedia Compression, Encryption, and Modeling](https://www.amazon.com/Cellular-Automata-Transforms-Applications-Compression/dp/0792378571)
* [Pseudorandom Pattern Generation by a 4-Neighborhood Cellular Automata Based on a Probabilistic Analysis](http://www.iaeng.org/publication/IMECS2008/IMECS2008_pp1908-1913.pdf)
* [A Study and Comparison of First and Second Order Cellular Automata with Examples](https://ntnuopen.ntnu.no/ntnu-xmlui/bitstream/handle/11250/258721/351877_FULLTEXT01.pdf?sequence=3)
* [Reversible CA](http://www.jcasim.de/main/node10.html)
