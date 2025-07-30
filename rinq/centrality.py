import networkx as nx
import numpy as np

def compute_eigenvector_centrality(G):
    return nx.eigenvector_centrality(G, max_iter=1000, tol=1e-6)

def compute_estrada_centrality(G):
    A = nx.to_numpy_array(G)
    evals, _ = np.linalg.eig(A)
    return np.sum(np.exp(evals)).real

def compute_estrada_centrality_manual(G):
    A = nx.to_numpy_array(G)
    n = A.shape[0]
    ec = np.zeros(n)
    for k in range(100):  # Series expansion
        term = np.linalg.matrix_power(A, k) / np.math.factorial(k)
        ec += np.diag(term)
    return ec
