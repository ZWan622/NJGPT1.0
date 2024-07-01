import numpy as np

# Define the function that removes the empty space
def remove_gaps(sequences):

    sample_seq = next(iter(sequences.values()))
    length = len(sample_seq)
    valid_positions = [i for i in range(length) if all(seq[i] != '-' and seq[i] != '?' for seq in sequences.values())]
    new_sequences = {key: ''.join(seq[i] for i in valid_positions) for key, seq in sequences.items()}
    return new_sequences, len(valid_positions)

# Define the function that transforms U to T
def convert_u_to_t(sequences):

    new_sequences = {key: seq.replace('U', 'T') for key, seq in sequences.items()}
    return new_sequences

# Define the function that calculates the Hamming distance
def calculate_hamming_distance(seq1, seq2):
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2)) / len(seq1)

# Define the function that constructs the distance matrix (Hamming distance)
def create_hamming_distance_matrix(sequences):
    seq_names = list(sequences.keys())
    num_seqs = len(seq_names)
    distance_matrix = np.zeros((num_seqs, num_seqs))

    # Remove empty space and transform U to T
    sequences = convert_u_to_t(sequences)
    sequences, valid_length = remove_gaps(sequences)

    for i in range(num_seqs):
        for j in range(i + 1, num_seqs):
            distance = calculate_hamming_distance(sequences[seq_names[i]], sequences[seq_names[j]])
            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance

    return distance_matrix, seq_names

def compute_ui(distance_matrix):
    n = len(distance_matrix)  # number of tips
    u = np.zeros(n)  # initialize u array

    for i in range(n):
        # Sum of distances from tip i to all other tips
        sum_distances = sum(distance_matrix[i][j] for j in range(n) if i != j)
        u[i] = sum_distances / (n - 2) if (n - 2) > 0 else 0

    return u

# Input sequence
sequences = {}


# Calculate the distance matrix (Hamming distance)
distance_matrix, seq_names = create_hamming_distance_matrix(sequences)

print("Distance Matrix (Hamming):")
print(distance_matrix)

u_values = compute_ui(distance_matrix)
print("u_i Values (Hamming):", u_values)
