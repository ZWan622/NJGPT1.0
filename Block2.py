##Block2
def compute_diff(distance_matrix, u_values):
    n = len(distance_matrix)  # number of tips
    diff_matrix = np.zeros((n, n))  # initialize diff matrix

    for i in range(n):
        for j in range(n):
            diff_matrix[i][j] = distance_matrix[i][j] - u_values[i] - u_values[j]

    return diff_matrix


# Compute DIFF values
diff_matrix = compute_diff(distance_matrix, u_values)

# Print DIFF matrix
print(diff_matrix)