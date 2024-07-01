import numpy as np
import math

# Define the function that removes the empty space
def remove_gaps(sequences):
    sample_seq = next(iter(sequences.values()))
    length = len(sample_seq)
    valid_positions = [i for i in range(length) if all(seq[i] != '-' and seq[i] != '?' for seq in sequences.values())]
    new_sequences = {key: ''.join(seq[i] for i in valid_positions) for key, seq in sequences.items()}
    return new_sequences, len(valid_positions)

# Define the function that transform U to T
def convert_u_to_t(sequences):
    new_sequences = {key: seq.replace('U', 'T') for key, seq in sequences.items()}
    return new_sequences

# Define the function that calculates the K2P distance
def compute_k2p_distances(seq1, seq2):
    transitions = 0
    transversions = 0
    for c1, c2 in zip(seq1, seq2):
        if c1 != c2:
            if (c1 in 'AG' and c2 in 'AG') or (c1 in 'CT' and c2 in 'CT'):
                transitions += 1
            else:
                transversions += 1
    p = transitions / len(seq1)
    q = transversions / len(seq1)

    w1 = 1 - 2 * p - q
    w2 = 1 - 2 * q

    try:
        d = -0.5 * math.log(w1) - 0.25 * math.log(w2)
        s = -0.5 * math.log(w1)
        v = -0.25 * math.log(w2)
        R = s / v if v != 0 else float('inf')
    except ValueError:
        d, s, v, R = float('inf'), float('inf'), float('inf'), float('inf')

    return d, s, v, R

# Define the function that constructs the distance matrix (K2p model)
def create_k2p_distance_matrix(sequences):
    seq_names = list(sequences.keys())
    num_seqs = len(seq_names)
    distance_matrix = np.zeros((num_seqs, num_seqs))

    sequences = convert_u_to_t(sequences)
    sequences, valid_length = remove_gaps(sequences)

    for i in range(num_seqs):
        for j in range(i + 1, num_seqs):
            d, s, v, R = compute_k2p_distances(sequences[seq_names[i]], sequences[seq_names[j]])
            distance_matrix[i, j] = d
            distance_matrix[j, i] = d

    return distance_matrix, seq_names


def compute_ui(distance_matrix):
    n = len(distance_matrix)
    u = np.zeros(n)

    for i in range(n):
        sum_distances = sum(distance_matrix[i][j] for j in range(n) if i != j)
        u[i] = sum_distances / (n - 2) if (n - 2) > 0 else 0

    return u

# Input sequence
sequences = {}

# Compute K2P distance matrix
distance_matrix, seq_names = create_k2p_distance_matrix(sequences)

print("Distance Matrix (K2P):")
print(distance_matrix)

u_values = compute_ui(distance_matrix)
print("u_i Values (K2P):", u_values)
