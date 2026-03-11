import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import KMeans
#from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix







def IdentifyNich(adata, group_by, cluster_name="niches", neighbors_k=30, niches_k=4):
    # Extract cell indices for adata
    cells = adata.obs.index
    group_labels = adata.obs.loc[cells, group_by].values

    # Identify unique groups
    groups = np.unique(group_labels)

    # Initialize the cell type matrix
    cell_type_mtx = np.zeros((len(cells), len(groups)), dtype=int)
    cell_idx_map = {cell: idx for idx, cell in enumerate(cells)}
    group_idx_map = {group: idx for idx, group in enumerate(groups)}

    for cell, group in zip(cells, group_labels):
        cell_idx = cell_idx_map[cell]
        group_idx = group_idx_map[group]
        cell_type_mtx[cell_idx, group_idx] = 1

    # Get spatial coordinates
    coords = adata.obs[['array_row','array_col']][adata.obs_names.isin(cells)]

    # Find nearest neighbors
    nbrs = NearestNeighbors(n_neighbors=neighbors_k).fit(coords)
    distances, indices = nbrs.kneighbors(coords)

    # Create sparse adjacency matrix
    adjacency = lil_matrix((len(cells), len(cells)))
    for i, neighbors in enumerate(indices):
        adjacency[i, neighbors] = 1
    adjacency = adjacency.tocsr()  # Convert to efficient CSR format

    # Calculate the sum matrix
    sum_mtx = adjacency @ cell_type_mtx

    # Scale data----- to center and scale each cell type across all individual cells
    scaled_data = (sum_mtx - sum_mtx.mean(axis=0)) / (sum_mtx.std(axis=0))

    # Perform k-means clustering
    kmeans = KMeans(n_clusters=niches_k, n_init=30, random_state=0)
    clusters = kmeans.fit_predict(scaled_data)

    # Store cluster labels in AnnData
    adata.obs[cluster_name] = pd.Categorical([clusters[cell_idx_map[cell]] for cell in adata.obs.index])

    return adata


#####How to use it:
#####sys.path.insert(0, './code/sourcePython/Spatial')  # adding marker Folder to the system path
#####import IdentifyCellNich2
#####adata=IdentifyCellNich2.IdentifyNich(adata=adata, xxxx)
#####



