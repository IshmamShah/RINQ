import numpy as np
import networkx as nx
from Bio import PDB
from Bio.PDB import PDBList

def fetch_pdb(pdb_id):
    pdbl = PDBList()
    pdb_file = pdbl.retrieve_pdb_file(pdb_id, file_format='pdb')
    return pdb_file

def calculate_distance(residue_i, residue_j):
    try:
        atom_i = residue_i["CA"].coord
        atom_j = residue_j["CA"].coord
        return np.linalg.norm(atom_i - atom_j)
    except KeyError:
        return float("inf")

def construct_protein_residue_network(pdb_id, interaction_cutoff=8.0):
    pdb_path = fetch_pdb(pdb_id)
    parser = PDB.PDBParser()
    structure = parser.get_structure("protein", pdb_path)
    G = nx.Graph()
    residues = []

    for chain in structure.get_chains():
        for residue in chain:
            if PDB.is_aa(residue):
                residue_id = residue.get_id()[1]
                residues.append((residue_id, residue))

    for i in range(len(residues)):
        for j in range(i + 1, len(residues)):
            dist = calculate_distance(residues[i][1], residues[j][1])
            if dist < interaction_cutoff:
                G.add_edge(residues[i][0], residues[j][0])

    return G

def create_adjacency_matrix(G, residue_order=None):
    if residue_order is None:
        residue_order = list(G.nodes())
    return nx.to_numpy_array(G, nodelist=residue_order)
