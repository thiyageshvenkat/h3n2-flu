# SeqKit Playground 🧬

### 1. Get Basic Statistics
```bash
seqkit stats results/split_sequences/*.fasta --all
```

### 2. View the Header/ID of Sequences
```bash
seqkit head -n 5 results/split_sequences/*.fasta | seqkit seq --name
```

### 3. Translate DNA to Protein (6-frame)
```bash
seqkit translate results/split_sequences/*.fasta > results/translated_proteins.faa
```

### 4. Search for a Motif (Regular Expression)
```bash
seqkit locate -p "ATG.{3,6}G" results/split_sequences/*.fasta
```

### 5. Sort Sequences by Length
```bash
seqkit sort --by-length --reverse results/split_sequences/*.fasta | seqkit head -n 10
```

### 6. Convert to Tabular Format
```bash
seqkit fx2tab results/split_sequences/*.fasta > results/sequences_table.tsv
```

### 7. Find Duplicates
```bash
seqkit rmdup --by-seq results/split_sequences/*.fasta --out-file results/clean_sequences.fasta
```

### 8. Find Open Reading Frames (ORFs)
```bash
seqkit translate results/split_sequences/*.fasta | seqkit seq --min-len 100 > results/proteins_only.faa
```

### 9. Rename sequences (Headers only)
```bash
seqkit replace -p "\|.*" -r "" results/split_sequences/*.fasta
```

### 10. Filter by GC Content (> 40%)
```bash
seqkit fx2tab --n --g results/split_sequences/*.fasta | awk '$2 > 40' | seqkit tab2fx
```

### 11. Mask poly-A tails with Ns
```bash
seqkit replace -i -p "A{10,}" -r "NNNNNNNNNN" results/split_sequences/*.fasta
```

---

> [!TIP]
> Most `seqkit` commands support **piping**.
> `seqkit grep -p "Spike" results/split_sequences/*.fasta | seqkit stats`