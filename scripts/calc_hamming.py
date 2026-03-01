from Bio import SeqIO
import sys

# hamming distance formula
def hamming_dist(s1, s2):
    if len(s1) != len(s2):
        print(f"Warning: Sequence lengths differ ({len(s1)} vs {len(s2)})", file=sys.stderr)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

# load sequences
ref = list(SeqIO.parse(sys.argv[1], "pdb-seqres"))[0]
var = list(SeqIO.parse(sys.argv[2], "fasta"))[0]

dist = hamming_dist(ref.seq, var.seq)

# We print only the number so Nextflow can easily capture it as a variable
print(dist)