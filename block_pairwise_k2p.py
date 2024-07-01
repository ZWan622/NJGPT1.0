import numpy as np
import math

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
    pairwise_sequences = remove_gaps_pairwise(sequences)

    for i in range(num_seqs):
        for j in range(i + 1, num_seqs):
            seq1, seq2 = pairwise_sequences[seq_names[i]][seq_names[j]]
            d, s, v, R = compute_k2p_distances(seq1, seq2)
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

# Input sequences
sequences = {}

# Compute distance matrix and u_i values
distance_matrix, seq_names = create_k2p_distance_matrix(sequences)
u_values = compute_ui(distance_matrix)


# Print the results
print("Distance Matrix (K2P):")
print(distance_matrix)

u_values_k2p = compute_ui(distance_matrix)
print("u_i Values (K2P):", u_values)

