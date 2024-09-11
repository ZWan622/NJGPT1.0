##Block5
def phylogenetic_reduction(distance_matrix, labels):
    while len(distance_matrix) > 2:
        # Step 1: Compute u_i for each node
        u_values = compute_ui(distance_matrix)

        # Step 2: Compute DIFF matrix
        diff_matrix = compute_diff(distance_matrix, u_values)

        # Step 3: Find the pair with the smallest DIFF and compute branch lengths
        min_pair = find_min_diff_pair(diff_matrix)
        i, j = min_pair

        # Compute branch lengths
        v_i = (distance_matrix[i][j] + u_values[i] - u_values[j]) / 2
        v_j = (distance_matrix[i][j] + u_values[j] - u_values[i]) / 2

        # Check for negative branch lengths and adjust
        if v_i < 0:
            v_j = v_j + v_i
            v_i = 0  # Set negative branch length to zero
        if v_j < 0:
            v_i = v_j + v_i
            v_j = 0  # Set negative branch length to zero

        # Step 4: Update the distance matrix
        distance_matrix = update_distance_matrix(distance_matrix, i, j)

        # Update labels
        new_label = f"({labels[i]}:{v_i:.4f},{labels[j]}:{v_j:.4f})"
        labels = [labels[k] for k in range(len(labels)) if k != i and k != j] + [new_label]

    # Final combination
    final_i, final_j = 0, 1
    final_vi = round((distance_matrix[final_i][final_j] + u_values[final_i] - u_values[final_j]) / 2, 4)
    final_vj = round((distance_matrix[final_i][final_j] + u_values[final_j] - u_values[final_i]) / 2, 4)

    # Modify final branch lengths to sum into one if one is non-zero
    if final_vi != 0 or final_vj != 0:
        final_vi += final_vj
        final_vj = 0

    final_label = f"({labels[final_i]}:{final_vi:.4f},{labels[final_j]}:{final_vj:.4f})"

    return final_label

# Initial labels for the nodes
labels = list(sequences.keys())

# Running the phylogenetic tree reduction
final_tree = phylogenetic_reduction(distance_matrix, labels)

print("final_tree:", final_tree)
