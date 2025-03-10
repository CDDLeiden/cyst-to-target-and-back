from typing import List, Optional, Union

import numpy as np
from sklearn.neighbors import NearestNeighbors


def compute_knn_distance(
    X: np.ndarray,
    Y: Optional[Union[np.ndarray, List[np.ndarray]]] = None,
    dist_metric: Optional[str] = None,
    k: int = 5,
    n_jobs: Optional[int] = None,
):
    """
    Computes the mean k-Nearest Neighbors distance between a vector of samples (X)
    and a set of query samples (y).

    Args:
        X: The set of samples that form kNN candidates
        Y: The samples for which to find the kNN for. If None, will find kNN in `X`
        dist_metric: The pairwise distance metric to define the neighborhood
        k: The number of neighbors to find
        n_jobs: Controls the parallelization
    """
    knn = NearestNeighbors(n_neighbors=k, metric=dist_metric, n_jobs=n_jobs)
    knn.fit(X)

    if not isinstance(Y, list):
        Y = [Y]

    distances, indices = [], []
    for queries in Y:
        if np.array_equal(X, queries):
            # Use k + 1 and filter out the first
            # because the sample will always be its own neighbor
            dist, ind = knn.kneighbors(queries, n_neighbors=k + 1)
            dist, ind = dist[:, 1:], ind[:, 1:]
        else:
            dist, ind = knn.kneighbors(queries, n_neighbors=k)

        distances.append(dist)
        indices.append(ind)

    if len(distances) == 1:
        assert len(indices) == 1
        distances = distances[0]
        indices = indices[0]

    return distances, indices
