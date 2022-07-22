# Progressive method
The multiple sequence alignment (MSA) is first performed with the progressive method with a bare bones implementation of the [Clustal](https://en.wikipedia.org/wiki/Clustal) algorithm.
The similarity matrix is given by the scores of pairwise Needleman-Wunsch alignments and an [UPGMA](https://en.wikipedia.org/wiki/UPGMA) tree is then built on top of that to get the order in which the sequences will be aligned, starting with the most similar pair and finishing with the most dissimilar sequence in the group.

## Multiple alignment class
The multiple alignments are stored in a vector of strings containing all the characters **at each position across sequences**, to reduce allocations and speed up iterations for consensus retrieval, score calculation and insertion at specific position.
So far I am not prioritizing  a base or amino acid when I have ties in the consensus calculation, the book does it in alphabetical order which seems very arbitrary.

