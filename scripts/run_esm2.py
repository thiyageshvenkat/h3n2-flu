import os
import sys
import importlib

os.environ.setdefault("HOME", "/tmp")
os.environ.setdefault("TORCHINDUCTOR_CACHE_DIR", "/tmp/torchinductor_cache")
os.environ.setdefault("TORCH_HOME", "/tmp/torch_home")
os.environ.setdefault("HF_HOME", "/tmp/huggingface")
os.environ.setdefault("TRANSFORMERS_CACHE", "/tmp/huggingface/hub")

import sentry_sdk

sentry_sdk.init(dsn=os.environ.get("SENTRY_DSN"), traces_sample_rate=0.0)

torch = importlib.import_module("torch")
SeqIO = importlib.import_module("Bio.SeqIO")
transformers = importlib.import_module("transformers")
AutoTokenizer = transformers.AutoTokenizer
AutoModelForMaskedLM = transformers.AutoModelForMaskedLM


def terminate_process(message, exit_code=1):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(exit_code)


"""
TODO: GPU, Mean-Pooling features
"""

# args
if len(sys.argv) != 2:
    terminate_process("Usage: python run_esm2.py <variant.fasta>")
fasta_input_path = sys.argv[1]
model_id = "facebook/esm2_t33_650M_UR50D"

# Load sequence
try:
    sequence_record = next(SeqIO.parse(fasta_input_path, "fasta"))
    protein_sequence = str(sequence_record.seq)
except StopIteration:
    terminate_process(f"No sequence found in {fasta_input_path}")

# Load model and tokenizer, tokenize, run inference
tokenizer = AutoTokenizer.from_pretrained(model_id)
model = AutoModelForMaskedLM.from_pretrained(model_id, output_hidden_states=True)
inputs = tokenizer(protein_sequence, return_tensors="pt")

with torch.no_grad():
    outputs = model(**inputs, labels=inputs["input_ids"])
    # Calculate perplexity
    perplexity_score = torch.exp(outputs.loss).item()
    # Extract latent embeddings [sequence_length, embedding_dim]
    embeddings_tensor = outputs.hidden_states[-1].squeeze(0)

# Define output filenames
perplexity_file = f"{fasta_input_path}.perplexity"
embeddings_file = f"{fasta_input_path}.embeddings.pt"

# Save outputs
with open(perplexity_file, "w") as f:
    f.write(str(perplexity_score))

torch.save(embeddings_tensor, embeddings_file)

print(embeddings_file)  # output for machine
print(
    f"Log: Processed {fasta_input_path}. Perplexity: {perplexity_score:.4f}",
    file=sys.stderr,
)  # output for user
