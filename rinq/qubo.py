import numpy as np
from scipy.optimize import dual_annealing


def calculate_d_hat(adj_matrix):
    """
    Compute normalized degree vector d_hat.

    Parameters:
        adj_matrix (np.ndarray): Adjacency matrix of the network.

    Returns:
        np.ndarray: Normalized degree vector.
    """
    d = np.sum(adj_matrix, axis=1)
    return d / np.linalg.norm(d)


def construct_qubo_matrix(adj_matrix, tau=5, P0=None, P1=None):
    """
    Construct a QUBO matrix to identify top-τ central nodes using Estrada-based centrality.

    Parameters:
        adj_matrix (np.ndarray): Adjacency matrix of the network.
        tau (int): Desired number of top nodes to select.
        P0 (float): Centrality weighting parameter. Defaults to 1/sqrt(n).
        P1 (float): Constraint penalty parameter. Defaults to 10 * n.

    Returns:
        np.ndarray: The QUBO matrix.
    """
    n = adj_matrix.shape[0]
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


def qubo_cost(x, Q):
    """
    Compute the QUBO objective function value for a given binary vector.

    Parameters:
        x (np.ndarray): Binary solution vector.
        Q (np.ndarray): QUBO matrix.

    Returns:
        float: Cost value.
    """
    return x.T @ Q @ x


def solve_qubo_simulated_annealing(Q, tau=5, num_reads=1000):
    """
    Solve the QUBO problem using classical simulated annealing.

    Parameters:
        Q (np.ndarray): QUBO matrix.
        tau (int): Desired number of selected nodes.
        num_reads (int): Number of independent annealing attempts.

    Returns:
        x_opt (np.ndarray): Binary solution with τ selected nodes.
        cost (float): Final energy of the best valid solution.
    """
    n = Q.shape[0]
    best_x = None
    best_cost = float("inf")

    def objective(x):
        x_bin = np.round(x).astype(int)
        if np.sum(x_bin) != tau:
            return 1e6  # Penalize invalid solutions
        return qubo_cost(x_bin, Q)

    bounds = [(0, 1)] * n

    for _ in range(num_reads):
        result = dual_annealing(objective, bounds, maxiter=100)
        x_bin = np.round(result.x).astype(int)
        if np.sum(x_bin) == tau:
            cost = qubo_cost(x_bin, Q)
            if cost < best_cost:
                best_cost = cost
                best_x = x_bin

    return best_x, best_cost
