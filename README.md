# RINQ – Residue Interaction Network Quantum engine

**RINQ** is a Python package that applies quantum and quantum-inspired algorithms to compute centrality measures in protein residue interaction networks.

## Installation

```bash
pip install .
```

## Usage Example

```python
from rinq.network import construct_protein_residue_network, create_adjacency_matrix
from rinq.centrality import compute_eigenvector_centrality
from rinq.qubo import construct_qubo_matrix, solve_qubo_simulated_annealing
```

# RINQ – Quantum Engine for Protein Network Centrality

**RINQ** is a quantum-enhanced method for identifying critical residues in protein networks, enabling scalable protein structure analysis for protein engineering and drug discovery.

---

## Contents

- **PDF files**: Contain the Jupyter Notebook code for the implementation of the RINQ pipeline for each of the proteins. All required packages are listed within these notebook files. A user can copy the notebook code and run it as-is to replicate the pipeline on their machine.
- **PNG files**: Contain the network visualizations generated during analysis.

---

## Reference

All theoretical details of the project are available in the following manuscript:

> [RinQ: Predicting central sites in proteins on current quantum computers](https://chemrxiv.org/engage/chemrxiv/article-details/686581403ba0887c333e255c).

---

## Contact

For questions or comments, please feel free to contact **smohtas@ncsu.edu** / **sishmam51@gmail.com**.

---


