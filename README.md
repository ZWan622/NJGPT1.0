# NJGPT1.0

##The relevant python code modules used by NJGPT, each with a different function

#block_complete_hamming.py:
Uses Hamming distance to measure the difference between sequences.
Applies complete deletion, meaning columns with missing data are removed entirely, ensuring every position has data.
Suitable for comparing highly similar DNA or protein sequences quickly, leading to the construction of a neighbor-joining (NJ) tree.

#block_complete_jc.py:
Uses the Jukes-Cantor (JC) model, a simple evolutionary model assuming equal substitution rates between all nucleotides.
Employs complete deletion, retaining only sites without any missing data.
This is useful for cases where sequences are relatively conserved, and it’s important to infer evolutionary distances based on the JC model.

#block_complete_k2p.py:
Applies the Kimura 2-Parameter (K2P) model, which accounts for different rates between transitions and transversions.
Also uses complete deletion, discarding any sites with missing data.
K2P is well-suited for situations where the evolutionary process includes unequal rates of transitions and transversions, providing a more accurate evolutionary distance.

#block_pairwise_hamming.py:
Similar to block_complete_hamming.py, it uses Hamming distance, but with pairwise deletion. In this method, only the positions where both sequences have data are used in the distance calculation.
This approach retains more data in cases with missing values, allowing for more comprehensive analysis while still constructing the NJ tree.

#block_pairwise_jc.py:
Uses the Jukes-Cantor model for evolutionary distance calculation, with pairwise deletion, meaning only sites where both sequences have data are considered.
This method is suitable for analyzing sequences with some missing data while maintaining accuracy in evolutionary distance estimation.

#block_pairwise_k2p.py:
Employs the Kimura 2-Parameter model, considering transitions and transversions, and uses pairwise deletion to handle missing data.
This model and method combination is useful when some sequences contain missing data, and the evolutionary process involves different substitution rates for transitions and transversions.

#Block2.py:
Purpose: Calculates the pairwise distance matrix between sequences.
Output: A distance matrix used to build the phylogenetic tree.

#Block3.py:
Purpose: Constructs a Neighbor-Joining (NJ) tree from the distance matrix.
Output: An initial tree structure showing evolutionary relationships.

#Block4.py:
Purpose: Refines the tree by adjusting branch lengths for better accuracy.
Output: A more accurate tree with optimized branch lengths.

#Block5.py:
Purpose: Exports the final tree in Newick format.
Output: A tree in Newick format, ready for visualization or further use.

#Block6.py:
Purpose: Visualizes the phylogenetic tree.
Output: A graphical representation of the tree.
