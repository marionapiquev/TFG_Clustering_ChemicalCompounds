# TFG: Characterization of commercial chemical compounds with pharmacological activity using clustering techniques.
In this repository, I provide the code used to carry out my study on "Characterization of commercial chemical compounds with pharmacological activity using clustering techniques". You will find a variety of clustering algorithms implemented in R. These algorithms include hierarchical clustering, k-means clustering, and density-based clustering. In addition, you will find scripts for data preprocessing and visualization, as well as datasets.

## Project Files Description

This project includes two directories, one for Scripts and another for Data.

#### Scripts:
- `get_molecular_descriptors.R`: Includes all functions required to extract molecular descriptors.
- `data_preprocessing.R`: Script with data transformation process implemented to prepare data for clustering.
- `KMeans_Clustering.R`: Includes a variety of functions to perform, analyze and visualize the clustering results of K-Means using molecular descriptors.
- `HC_Clustering.R`: Includes a variety of functions to perform, analyze and visualize the clustering results of Hierarchical Clustering using molecular descriptors.
- `clustering_visualitzation_dbscan.R`: Includes a variety of functions to analyze and visualize the clustering results of DBSCAN.
- `clustering_visualitzation_hdbscan.R`: Includes a variety of functions to analyze and visualize the clustering results of HDBSCAN.
- `clustering_visualitzation_optics.R`: Includes a variety of functions to analyze and visualize the clustering results of OPTICS.
- `Clustering_Fingerprints.R`: Includes a variety of clustering algorithms implemented to cluster chemical compounds based on their fingerprints. These algorithms are hierarchical clustering, K-means, MDS+HC and MDS+K-Means.
- `Validation.ipynb`: Jupiter Notebook that includes the supervised learning validation for question 1.
#### Data:
- `df_molecular_descriptors.csv`: Dataframe obtained using *extract_molecular_descriptors* function from `get_molecular_descriptors.R` file. It contains the molecular descriptors of 9213 chemical compounds.
- `df_molecular_descriptors_scaled.csv`: Dataframe obtained executing the `data_preprocessing.R` file. It contains the molecular descriptors of 9213 chemical compounds after data preprocessing.
- `df_smiles.csv`: Dataframe obtained using *extract_smiles* function from `get_molecular_descriptors.R` file. It contains the smiles of 9213 chemical compounds.
- `df_info_mols.csv`: Dataframe that contains additional information about each chemical compound such as generic name and phase.

The original SDF file can't be uploaded due to licensing reasons.
