import sys
from Bio import SeqIO

def calc_cpg(fasta_file):
    # Load sequence
    record = next(SeqIO.parse(fasta_file, "fasta"))
    seq = str(record.seq).upper()
    if len_seq == 0 or count_c == 0 or count_g == 0:
        cpg_ratio = 0.0 # extra safety even though already should be checked
    # Calculate base counts
    len_seq = len(seq)
    count_c = seq.count('C')
    count_g = seq.count('G')
    count_cg = seq.count('CG')
    
    # Calculate Observed vs Expected CpG Ratio
    if count_c == 0 or count_g == 0:
        cpg_ratio = 0.0
    else:
        expected_cg = (count_c * count_g) / len_seq
        cpg_ratio = count_cg / expected_cg
        
    # Output
    print(cpg_ratio) # Output for machine
    print(f"Log: CpG Ratio: {cpg_ratio}", file=sys.stderr) # Output for user

if __name__ == "__main__":
    # python calc_cpg.py variant.fasta > cpg_score.txt
    calc_cpg(sys.argv[1])