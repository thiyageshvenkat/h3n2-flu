import os
import sys
import sentry_sdk
import numpy as np
from Bio import SeqIO
from sklearn.metrics import mutual_info_score
from itertools import combinations

sentry_sdk.init(dsn=os.environ.get("SENTRY_DSN"), traces_sample_rate=0.0)



def calculate_epistasis_matrix(fasta_file, threshold=0.05):
    # Load aligned sequences into 2D matrix
    sequence_records = list(SeqIO.parse(fasta_file, "fasta"))
    max_sequence_length = max(len(record.seq) for record in sequence_records)
    sequence_matrix = np.array([list(str(record.seq).ljust(max_sequence_length, "-")) for record in sequence_records])

    # Find mutation sites
    variant_positions = [
        column_index
        for column_index in range(sequence_matrix.shape[1])
        if len(set(sequence_matrix[:, column_index]) - {"-", "X"}) > 1
    ]

    # CSV header
    print("Site1,Site2,Mutual_Information")

    # Calculate Mutual Information for every pair of mutations
    for position_1, position_2 in combinations(variant_positions, 2):
        mutual_info_value = mutual_info_score(sequence_matrix[:, position_1], sequence_matrix[:, position_2])
        # Only output significant epistatic links
        if mutual_info_value > threshold:
            print(f"{position_1},{position_2},{mutual_info_value:.4f}")


if __name__ == "__main__":
    # python calc_epistasis.py aligned_bundle.fasta > epistasis_matrix.csv
    fasta_input_path = sys.argv[1]
    score_threshold = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
    calculate_epistasis_matrix(fasta_input_path, score_threshold)