from Bio import SeqIO
import sys

# hamming distance formula
def calculate_hamming_distance(sequence_1, sequence_2):
    if len(sequence_1) != len(sequence_2):
        print(
            f"Warning: Sequence lengths differ ({len(sequence_1)} vs {len(sequence_2)})",
            file=sys.stderr,
        )
    return sum(char_1 != char_2 for char_1, char_2 in zip(sequence_1, sequence_2))


# load sequences
reference_record = next(SeqIO.parse(sys.argv[1], "pdb-seqres"))
variant_record = next(SeqIO.parse(sys.argv[2], "fasta"))
total_distance = calculate_hamming_distance(str(reference_record.seq), str(variant_record.seq))
print(total_distance) # Output for machine
print(f"Log: Calculated distance for {sys.argv[2]}", file=sys.stderr) # Output for user