##Block4

def update_distance_matrix(distance_matrix, i, j):
    n = len(distance_matrix)
    new_distances = []

    # Calculate distances from new node M to all other nodes k
    for k in range(n):
        if k != i and k != j:
            distance_to_k = (distance_matrix[i][k] + distance_matrix[j][k] - distance_matrix[i][j]) / 2
            new_distances.append(distance_to_k)

    # Create a new reduced matrix with new node
    new_matrix = []
    for k in range(n):
        if k != i and k != j:
            new_row = [distance_matrix[k][m] for m in range(n) if m != i and m != j]
            new_row.append(new_distances.pop(0))
            new_matrix.append(new_row)

    # Append the new node's distances to all other nodes
    new_matrix.append([row[-1] for row in new_matrix] + [0])  # The last entry is distance to itself (0)

    return new_matrix


# Update the distance matrix by replacing tips i and j with the new node
updated_matrix = update_distance_matrix(distance_matrix, i, j)
print(updated_matrix)