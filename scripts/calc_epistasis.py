import sys
import numpy as np
from Bio import SeqIO
from sklearn.metrics import mutual_info_score
from itertools import combinations

def calc_epistasis(fasta_file, threshold=0.05):
    # Load aligned sequences into 2D matrix
    records = list(SeqIO.parse(fasta_file, "fasta"))
    max_len = max(len(r.seq) for r in records)
    seqs = np.array([list(str(r.seq).ljust(max_len, '-')) for r in records])

    # Find mutation sites
    mut_sites = [i for i in range(seqs.shape[1]) if len(set(seqs[:, i]) - {'-', 'X'}) > 1]

    # Print CSV header
    print("Site1,Site2,Mutual_Information")

    # Calculate Mutual Information for every pair of mutations
    for site1, site2 in combinations(mut_sites, 2):
        mi = mutual_info_score(seqs[:, site1], seqs[:, site2])
        
        # Only output significant epistatic links
        if mi > threshold: 
            print(f"{site1},{site2},{mi:.4f}")

if __name__ == "__main__":
    # python calc_epistasis.py aligned_bundle.fasta > epistasis_matrix.csv
    fasta_in = sys.argv[1]
    thresh = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
    calc_epistasis(fasta_in, thresh)