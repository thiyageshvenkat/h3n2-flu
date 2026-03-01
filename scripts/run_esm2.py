import torch
import sys
from Bio import SeqIO
from transformers import AutoTokenizer, AutoModelForMaskedLM

# Load sequence
record = list(SeqIO.parse(sys.argv[1], "fasta"))[0]
sequence = str(record.seq)

# Load model
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModelForMaskedLM.from_pretrained(model_name, output_hidden_states=True)

inputs = tokenizer(sequence, return_tensors="pt")
with torch.no_grad():
    outputs = model(**inputs, labels=inputs["input_ids"])
    
    # Calculate perplexity
    perplexity = torch.exp(outputs.loss).item()
    
    # Extract latent embeddings
    # [1, sequence_length, embedding_dim]
    embeddings = outputs.hidden_states[-1].squeeze(0) 

# Save outputs for main.nf to get
with open(f"{sys.argv[1]}.perplexity", "w") as f:
    f.write(str(perplexity))

torch.save(embeddings, f"{sys.argv[1]}.embeddings.pt")