# Progressive method
The multiple sequence alignment (MSA) is first performed with the progressive method with a bare bones implementation of the [Clustal](https://en.wikipedia.org/wiki/Clustal) algorithm.
The similarity matrix is given by the scores of pairwise Needleman-Wunsch alignments and an [UPGMA](https://en.wikipedia.org/wiki/UPGMA) tree is then built on top of that to get the order in which the sequences will be aligned, starting with the most similar pair and finishing with the most dissimilar sequence in the group.


