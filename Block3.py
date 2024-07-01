##Block3
def find_min_diff_pair(diff_matrix):
    n = len(diff_matrix)
    min_value = float('inf')
    min_index = (-1, -1)

    for i in range(n):
        for j in range(i + 1, n):  # Consider only upper triangle for undirected graph
            if diff_matrix[i][j] < min_value:
                min_value = diff_matrix[i][j]
                min_index = (i, j)

    return min_index


# Find the pair with the smallest DIFF
min_pair = find_min_diff_pair(diff_matrix)
i, j = min_pair

# Calculate branch lengths
v_i = round((distance_matrix[i][j] + u_values[i] - u_values[j]) / 2, 4)
v_j = round((distance_matrix[i][j] + u_values[j] - u_values[i]) / 2, 4)

# Output results
print("Minimum DIFF pair:", min_pair)
print("Branch length v_i:", v_i)
print("Branch length v_j:", v_j)