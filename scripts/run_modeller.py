import sys
import os
from modeller import *
from modeller.automodel import *

# sys.argv[1] = variant sequence, sys.argv[2] = template pdb
seq_file, tpl_file = sys.argv[1], sys.argv[2]
tpl_name = os.path.splitext(os.path.basename(tpl_file))[0]

env = environ()
env.io.atom_files_directory = ['.']

# Alignment
aln = alignment(env)
aln.append_model(model(env, file=tpl_file), align_codes=tpl_name, atom_files=tpl_file)
aln.append(file=seq_file, align_codes='variant')
aln.align2d()
aln.write(file='alignment.ali', alignment_format='PIR')

# Build 3D model
a = automodel(env, alnfile='alignment.ali', knowns=tpl_name, sequence='variant')
a.starting_model = a.ending_model = 1
a.make()

# Rename safely using Modeller's exact output file
out_name = f"{os.path.splitext(os.path.basename(seq_file))[0]}_model.pdb"
os.rename(a.outputs[0]['name'], out_name)