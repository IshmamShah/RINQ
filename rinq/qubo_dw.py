import numpy as np
from dimod import BinaryQuadraticModel, SimulatedAnnealingSampler


def calculate_d_hat(adj_matrix):
    """Compute normalized degree vector d_hat."""
    d = np.sum(adj_matrix, axis=1)
    d_hat = d / np.linalg.norm(d)
    return d_hat


def construct_qubo_matrix_dw(adj_matrix, tau, P0=None, P1=None):
    """
    Construct the QUBO matrix using the formulation based on eigencentrality
    and Estrada centrality, adapted for quantum annealing.

    Parameters:
        adj_matrix (np.ndarray): The adjacency matrix.
        tau (int): Number of top nodes to select.
        P0 (float): Scaling parameter for the centrality term.
        P1 (float): Penalty parameter.

    Returns:
        Q (np.ndarray): QUBO matrix.
    """
    n = len(adj_matrix)
    if P0 is None:
        P0 = 1 / np.sqrt(n)
    if P1 is None:
        P1 = 10 * n

    d_hat = calculate_d_hat(adj_matrix)
    I = np.eye(n)
    U = np.ones((n, n)) - I
    C = (1 - 2 * tau) * I + U

    Q = -P0 * (adj_matrix @ np.outer(d_hat, d_hat) @ adj_matrix)
    Q += -P0 * (adj_matrix @ np.outer(d_hat, d_hat) @ adj_matrix @ adj_matrix)
    Q += P1 * C
    return Q


def solve_qubo_with_dwave(Q, tau, num_reads=10000, beta_range=(0.1, 4.0)):
    """
    Solve the QUBO using simulated annealing (D-Wave classical sampler).
    Filters solutions with exactly tau selected nodes.

    Parameters:
        Q (np.ndarray): QUBO matrix.
        tau (int): Desired number of selected nodes.
        num_reads (int): Number of annealing reads.
        beta_range (tuple): Annealing schedule range.

    Returns:
        top_nodes (List[int]): Node indices with value 1 in best solution.
        best_sample (dict): Raw binary sample.
        best_energy (float): Corresponding energy.
    """
    n = Q.shape[0]
    bqm = BinaryQuadraticModel.from_numpy_matrix(Q)
    sampler = SimulatedAnnealingSampler()
    response = sampler.sample(bqm, num_reads=num_reads, beta_range=beta_range)

    valid_solutions = []
    for sample, energy in response.data(['sample', 'energy']):
        if sum(sample.values()) == tau:
            valid_solutions.append((sample, energy))

    if not valid_solutions:
        return None, None, None

    best_sample, best_energy = min(valid_solutions, key=lambda x: x[1])
    top_nodes = [node for node, val in best_sample.items() if val == 1]

    return top_nodes, best_sample, best_energy
