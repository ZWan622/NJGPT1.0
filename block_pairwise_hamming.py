import numpy as np


# Define the function that removes the empty space
def remove_gaps_pairwise(sequences):
    sample_seq = next(iter(sequences.values()))
    length = len(sample_seq)
    pairwise_sequences = {}

    for key1, seq1 in sequences.items():
        pairwise_sequences[key1] = {}
        for key2, seq2 in sequences.items():
            if key1 != key2:
                valid_positions = [i for i in range(length) if seq1[i] != '-' and seq1[i] != '?' and seq2[i] != '-' and seq2[i] != '?']
                new_seq1 = ''.join(seq1[i] for i in valid_positions)
                new_seq2 = ''.join(seq2[i] for i in valid_positions)
                pairwise_sequences[key1][key2] = (new_seq1, new_seq2)

    return pairwise_sequences

# Define the function that constructs the distance matrix (Hamming distance)
def compute_distance_matrix(pairwise_sequences):
    keys = list(pairwise_sequences.keys())
    size = len(keys)
    distance_matrix = [[0] * size for _ in range(size)]

    key_index = {key: idx for idx, key in enumerate(keys)}

    for key1, pairs in pairwise_sequences.items():
        for key2, (seq1, seq2) in pairs.items():
            length = len(seq1)
            if length == 0:
                distance = 0
            else:
                different_nucleotide_count = sum(1 for a, b in zip(seq1, seq2) if a != b)
                distance = different_nucleotide_count / length
            idx1 = key_index[key1]
            idx2 = key_index[key2]
            distance_matrix[idx1][idx2] = distance

    return distance_matrix

# Define the function that transforms U to T
def convert_u_to_t(sequences):
    return {key: seq.replace('U', 'T') for key, seq in sequences.items()}

def compute_ui(distance_matrix):
    n = len(distance_matrix)  # number of tips
    u = np.zeros(n)  # initialize u array

    for i in range(n):
        # Sum of distances from tip i to all other tips
        sum_distances = sum(distance_matrix[i][j] for j in range(n) if i != j)
        u[i] = sum_distances / (n - 2) if (n - 2) > 0 else 0

    return u

# Example usage
sequences = {}

# Convert U to T
sequences = convert_u_to_t(sequences)


# Remove gaps pairwise
pairwise_sequences = remove_gaps_pairwise(sequences)

# Compute distance matrix
distance_matrix = compute_distance_matrix(pairwise_sequences)

print("Pairwise sequences after removing gaps:", pairwise_sequences)
print("Distance matrix:", distance_matrix)
u_values = compute_ui(distance_matrix)
print("u_i Values (Hamming):", u_values)

