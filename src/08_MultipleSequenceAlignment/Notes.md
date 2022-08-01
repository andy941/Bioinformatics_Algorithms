# Progressive method
The multiple sequence alignment (MSA) is first performed with the progressive method with a bare bones implementation of the [Clustal](https://en.wikipedia.org/wiki/Clustal) algorithm.
The similarity matrix is given by the scores of pairwise Needleman-Wunsch alignments and an [UPGMA](https://en.wikipedia.org/wiki/UPGMA) tree is then built on top of that to get the order in which the sequences will be aligned, starting with the most similar pair and finishing with the most dissimilar sequence in the group.
I am implementing this only for proteins for now.

## similarity
The similarity is calculated by pure pairwise alignment *score* for the purpose of this exercise even though it's far from optimal. Next chapter delves deeper into phylogenetic analysis. I am using 0,-1 and -1 as match, mismatch and gap scores as a rough estimate of sequence distance (-score()).

## Multiple alignment class
The multiple alignments are stored in a vector of strings containing all the characters **at each position across sequences**, to reduce allocations and speed up iterations for consensus retrieval, score calculation and insertion at specific position.
So far I am not prioritizing  a base or amino acid when I have ties in the consensus calculation, the book does it in alphabetical order which seems very arbitrary.

## UPGMA
To calculate the UPGMA efficiently two similarity matrices are created of fixed size, sim_mat and result. Each iteration the pointers to these matrices are swapped and the results matrix is zeroed  and used again for the next clustering step. Each iteration I'll perform a global alignment and then get the consensus from the accumulating multiple alignment.
