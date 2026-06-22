import glob
import os
import sys

import pandas as pd
import sentry_sdk

sentry_sdk.init(dsn=os.environ.get("SENTRY_DSN"), traces_sample_rate=0.0)


def terminate_process(error_message, exit_code=1):
    print(f"FATAL: {error_message}", file=sys.stderr)
    sys.exit(exit_code)


# Accept one wrapper directory or multiple feature directories.
if len(sys.argv) < 3:
    terminate_process(
        "Usage: python aggregate.py <input_directory> [<input_directory> ...] "
        "<output_dataset.parquet>"
    )

# The final argument is always the output; every preceding argument is an input folder.
input_directories = sys.argv[1:-1]
final_output_parquet = sys.argv[-1]

# CHECK: Verify that every input directory exists.
for input_directory in input_directories:
    if not os.path.isdir(input_directory):
        terminate_process(f"Input directory not found: {input_directory}")

# Dictionary to store all features
variant_data = {}


def get_or_create_variant_record(variant_identifier):
    if variant_identifier not in variant_data:
        variant_data[variant_identifier] = {}
    return variant_data[variant_identifier]


def extract_stability_record(filepath, line):
    if not line:
        print(f"WARNING: Skipping empty stability file {filepath}", file=sys.stderr)
        return None

    parts = line.split()
    if len(parts) < 2:
        print(
            f"WARNING: Skipping malformed stability file {filepath}: {line!r}",
            file=sys.stderr,
        )
        return None

    pdb_name = os.path.basename(parts[0])
    variant_identifier = pdb_name.removesuffix("__Repair.pdb").removesuffix(
        "_Repair.pdb"
    )
    return variant_identifier, parts[1]


# Function to load scalar features from files with redundancy checks and error handling
def load_scalar_features(
    results_directory, pattern, suffix, column_name, record_extractor=None
):
    for filepath in glob.glob(
        os.path.join(results_directory, "**", pattern), recursive=True
    ):
        filename = os.path.basename(filepath)

        with open(filepath, "r", encoding="utf-8") as file:
            line = file.read().strip()

        if record_extractor is None:
            variant_identifier = filename.removesuffix(suffix)
            raw_value = line
        else:
            extracted_record = record_extractor(filepath, line)
            if extracted_record is None:
                continue
            variant_identifier, raw_value = extracted_record

        variant_identifier = variant_identifier.removesuffix("_.split")

        try:
            value = float(raw_value)
        except ValueError:
            print(
                f"WARNING: Skipping invalid {column_name} in {filepath}: {raw_value!r}",
                file=sys.stderr,
            )
            continue

        variant_record = get_or_create_variant_record(variant_identifier)
        if column_name in variant_record:
            print(
                f"WARNING: Duplicate {column_name} for {variant_identifier}; "
                f"overwriting {variant_record[column_name]!r} with {value!r} "
                f"from {filepath}",
                file=sys.stderr,
            )

        variant_record[column_name] = value


print(
    f"Log: Scanning input directories: {', '.join(input_directories)}",
    file=sys.stderr,
)

for input_directory in input_directories:
    # Parse Hamming distances
    load_scalar_features(input_directory, "*_dist.txt", "_dist.txt", "Hamming_Distance")

    # Parse CpG ratios
    load_scalar_features(input_directory, "*_cpg.txt", "_cpg.txt", "CpG_Ratio")

    # Parse ESM-2 perplexity scores
    load_scalar_features(
        input_directory, "*.perplexity", ".fa.perplexity", "ESM2_Perplexity"
    )

    # Parse FoldX stability scores
    load_scalar_features(
        input_directory,
        "*.fxout",
        None,
        "Stability_Score",
        record_extractor=extract_stability_record,
    )

# CHECK: Ensure features were actually found
if not variant_data:
    terminate_process(f"No feature files found in: {', '.join(input_directories)}")

# Convert dictionary to pandas DataFrame
master_dataset = pd.DataFrame.from_dict(variant_data, orient="index")
master_dataset.index.name = "Variant_ID"
master_dataset = master_dataset.fillna(0.0)

# Save final output as .parquet
master_dataset.to_parquet(final_output_parquet, engine="pyarrow")

# final outputs
print(final_output_parquet)  # output for machine
print(
    f"Log: Aggregated {len(master_dataset)} variants into {final_output_parquet}",
    file=sys.stderr,
)  # output for user
