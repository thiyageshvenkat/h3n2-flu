import os
import sys
import sentry_sdk
import torch
import Bio.PDB
import numpy as np
from torch_geometric.data import Data

sentry_sdk.init(dsn=os.environ.get("SENTRY_DSN"), traces_sample_rate=0.0)


def build_protein_graph(pdb_path, embedding_path, output_path):
    # Load 3D structure
    parser = Bio.PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)

    coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:  # get Alpha Carbon for node pos
                    coords.append(residue["CA"].get_coord())

    node_coords = torch.tensor(np.array(coords), dtype=torch.float)

    # Load ESM-2 Embeddings
    # .pt file from RUN_ESM2
    node_features = torch.load(embedding_path)

    # create edges based on 3D distance
    # 8.0 Angstroms is the cutoff for residue contacts
    dist_matrix = torch.cdist(node_coords, node_coords)
    edge_index = (dist_matrix < 8.0).nonzero(as_tuple=False).t()

    # GNN Package
    protein_graph = Data(x=node_features, edge_index=edge_index, pos=node_coords)

    torch.save(protein_graph, output_path)
    print(f"Graph generated: {output_path}")


if __name__ == "__main__":
    build_protein_graph(sys.argv[1], sys.argv[2], sys.argv[3])
