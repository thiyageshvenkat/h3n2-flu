import os
import sys
import shutil
import glob
from modeller import *
from modeller.automodel import *

def terminate_process(error_message, exit_code=2):
    """Prints a clear error message to stderr and terminates the script."""
    print(f"ERROR: {error_message}", file=sys.stderr)
    sys.exit(exit_code)

# Ensure script has the correct number of arguments
if len(sys.argv) != 3:
    terminate_process("Usage: python build_model.py <variant_sequence.fasta> <template_structure.pdb>")

# Define alternative input paths
sequence_input_path = sys.argv[1]
template_pdb_path = sys.argv[2]

# CHECK: Verify that input files exist before starting the engine
if not os.path.isfile(sequence_input_path):
    terminate_process(f"Sequence file not found: {sequence_input_path}")
if not os.path.isfile(template_pdb_path):
    terminate_process(f"Template PDB not found: {template_pdb_path}")

# Extract directory and naming info to help Modeller find its files
template_directory = os.path.abspath(os.path.dirname(template_pdb_path) or ".")
template_filename = os.path.basename(template_pdb_path)
template_structure_id = os.path.splitext(template_filename)[0]
variant_identifier = os.path.splitext(os.path.basename(sequence_input_path))[0]

# Initialize the Modeller environment
log.verbose()
modeller_environment = environ()
# Ensures Modeller searches the template directory for the .pdb file
modeller_environment.io.atom_files_directory = ['.', template_directory]

# Align 4O5I.pdb to the viral sequence
try:
    reference_structure = model(modeller_environment, file=template_filename)
    sequence_alignment = alignment(modeller_environment)
    # Map the template structure to the alignment object
    sequence_alignment.append_model(reference_structure, align_codes=template_structure_id, atom_files=template_filename)
    # Add the target viral variant sequence
    sequence_alignment.append(file=sequence_input_path, align_codes='variant', alignment_format='FASTA')
    # Perform 2D alignment and generate PIR file for modeling
    sequence_alignment.align2d()
    sequence_alignment.write(file='alignment.ali', alignment_format='PIR')
    sequence_alignment.check()
    print(f"Log: Alignment completed for {variant_identifier}", file=sys.stderr)

except Exception as alignment_error:
    terminate_process(f"Alignment phase failed: {alignment_error}")

# Build the 3D model w/ comparative modeling
homology_modeler = automodel(
    modeller_environment,
    alnfile='alignment.ali',
    knowns=template_structure_id,
    sequence='variant',
    assess_methods=(assess.DOPE, assess.GA341)
)
homology_modeler.starting_model = 1
homology_modeler.ending_model = 1
homology_modeler.make()

# Identify and save the best model
successful_models = [model_info for model_info in homology_modeler.outputs if model_info.get('failure') is None]

if not successful_models:
    terminate_process("Modeller failed to produce a valid 3D structure. Check logs for geometry issues.")

# Select the model with the lowest (best) DOPE score
best_model_data = min(successful_models, key=lambda m: m.get('DOPE score', float('inf')))
generated_pdb_path = best_model_data['name']

# save final output
final_variant_pdb = f"{variant_identifier}_model.pdb"
shutil.copy(generated_pdb_path, final_variant_pdb)

# remove temp files
temporary_file_patterns = ['*.log', '*.ini', '*.rsr', '*.sch', '*.V9999*', 'variant.*', 'alignment.pap']
for pattern in temporary_file_patterns:
    for temporary_file in glob.glob(pattern):
        try:
            os.remove(temporary_file)
        except OSError:
            pass # Ignore files that are currently locked or already deleted

# final outputs
print(final_variant_pdb) # output for machine
print(f"Log: {variant_identifier} modeled. DOPE Score: {best_model_data.get('DOPE score')}", file=sys.stderr) # output for user