from sklearn.cluster import KMeans
import pandas as pd
import matplotlib.pyplot as plt


def __pop_coords(sequences, filtered_rows, n_clusters=2):
    """
    Create populations according to their coordinates (x,y) using KMeans clustering.
    
    Parameters:
    sequences (dict): A dictionary containing sequences keyed by index.
    filtered_rows (pd.DataFrame): A DataFrame containing 'x' and 'y' coordinate columns.
    n_clusters (int): Number of clusters for KMeans. Default is 2.
    
    Returns:
    tuple: Two lists representing the populations split by clusters.
    """

    # Create two empty populations
    populations = [[] for _ in range(n_clusters)]

    extracted_values = filtered_rows.iloc[:, :2].astype(float)
    coordinates = extracted_values[['x', 'y']].values

    # Perform KMeans clustering
    kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(coordinates)
    labels = kmeans.labels_

    # Iterate through extracted values and assign to populations based on clusters
    for index, label in zip(extracted_values.index, labels):
        sequence = sequences.get(str(index), None)  # Convert index to string to match sequence keys
        if sequence is not None:
            populations[label].append(sequence)
    # Plotting the clusters
    plt.figure(figsize=(10, 8))
    plt.scatter(coordinates[:, 0], coordinates[:, 1], c=labels, cmap='viridis', marker='o')
    centers = kmeans.cluster_centers_
    plt.scatter(centers[:, 0], centers[:, 1], c='red', s=300, alpha=0.6, marker='x')
    plt.title('KMeans Clustering of Coordinates')
    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.show()

    return populations

# Example usage
sequences = {
    '0': 'ATGC',
    '1': 'GCTA',
    '2': 'TACG',
    '3': 'CGTA'
}

filtered_rows = pd.DataFrame({
    'x': [1.0, 2.0, 3.0, 4.0],
    'y': [1.0, 2.0, 3.0, 4.0]
})

population_1, population_2 = __pop_coords(sequences, filtered_rows)
print("Population 1:", population_1)
print("Population 2:", population_2)
