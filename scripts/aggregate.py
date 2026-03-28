import os
import sys
import glob
import pandas as pd
import torch
import sentry_sdk

sentry_sdk.init(dsn=os.environ.get("SENTRY_DSN"), traces_sample_rate=0.0)


def terminate_process(error_message, exit_code=1):
    print(f"FATAL: {error_message}", file=sys.stderr)
    sys.exit(exit_code)

# Ensure script has the correct number of arguments
if len(sys.argv) != 3:
    terminate_process("Usage: python aggregate.py <results_directory> <output_dataset.parquet>")

# Define better input and output path name
results_directory = sys.argv[1]
final_output_parquet = sys.argv[2]

# CHECK: Verify that the results directory exists
if not os.path.isdir(results_directory):
    terminate_process(f"Results directory not found: {results_directory}")

# Dictionary to store all features
variant_data = {}

def get_or_create_variant_record(variant_identifier):
    if variant_identifier not in variant_data:
        variant_data[variant_identifier] = {}
    return variant_data[variant_identifier]

print(f"Log: Scanning {results_directory} for feature files...", file=sys.stderr)

# Parse Hamming distances
for filepath in glob.glob(os.path.join(results_directory, "**", "*__dist.txt"), recursive=True):
    filename = os.path.basename(filepath)
    variant_identifier = filename.replace("__dist.txt", "")
    with open(filepath, "r") as file:
        get_or_create_variant_record(variant_identifier)["Hamming_Distance"] = float(file.read().strip())

# Parse CpG ratios
for filepath in glob.glob(os.path.join(results_directory, "**", "*_cpg.txt"), recursive=True):
    filename = os.path.basename(filepath)
    variant_identifier = filename.replace("_cpg.txt", "")
    with open(filepath, "r") as file:
        get_or_create_variant_record(variant_identifier)["CpG_Ratio"] = float(file.read().strip())

# Parse ESM-2 perplexity scores
for filepath in glob.glob(os.path.join(results_directory, "**", "*.perplexity"), recursive=True):
    filename = os.path.basename(filepath)
    variant_identifier = filename.replace(".perplexity", "")
    with open(filepath, "r") as file:
        get_or_create_variant_record(variant_identifier)["ESM2_Perplexity"] = float(file.read().strip())

# Parse ESM-2 embeddings
for filepath in glob.glob(os.path.join(results_directory, "**", "*.embeddings.pt"), recursive=True):
    filename = os.path.basename(filepath)
    variant_identifier = filename.replace(".embeddings.pt", "")
    embedding_tensor = torch.load(filepath, weights_only=True)
    get_or_create_variant_record(variant_identifier)["ESM2_Embeddings"] = embedding_tensor.tolist()

# CHECK: Ensure features were actually found
if not variant_data:
    terminate_process(f"No feature files found in {results_directory}.")

# Convert dictionary to pandas DataFrame
master_dataset = pd.DataFrame.from_dict(variant_data, orient="index")
master_dataset.index.name = "Variant_ID"
master_dataset = master_dataset.fillna(0.0)

# Save final output as .parquet
master_dataset.to_parquet(final_output_parquet, engine="pyarrow")

# final outputs
print(final_output_parquet) # output for machine
print(f"Log: Aggregated {len(master_dataset)} variants into {final_output_parquet}", file=sys.stderr) # output for user