import sys
from Bio import SeqIO

def calc_cpg(fasta_file):
    # Load sequence
    record = list(SeqIO.parse(fasta_file, "fasta"))[0]
    seq = str(record.seq).upper()
    
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
    print(cpg_ratio)

if __name__ == "__main__":
    # python calc_cpg.py variant.fasta > cpg_score.txt
    calc_cpg(sys.argv[1])