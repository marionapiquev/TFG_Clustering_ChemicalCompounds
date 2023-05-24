# TFG: Characterization of commercial chemical compounds with pharmacological activity using clustering techniques.
In this repository, I provide the code used to carry out my study on "Characterization of commercial chemical compounds with pharmacological activity using clustering techniques". You will find a variety of clustering algorithms implemented in R. These algorithms include hierarchical clustering, k-means clustering, and density-based clustering. In addition, you will find scripts for data preprocessing and visualization, as well as datasets.

## Project Files Description

This project includes two directories, one for Scripts and another for Data.

#### Scripts:
- `get_molecular_descriptors.R`: Includes all functions required to extract molecular descriptors.
- `data_preprocessing.R`: Script with data transformation process implemented to prepare data for clustering.

#### Data:
- `df_molecular_descriptors.csv`: Dataframe obtained using *extract_molecular_descriptors* function from `get_molecular_descriptors.R` file. It contains the molecular descriptors of 9213 chemical compounds.
- `df_molecular_descriptors_scaled.csv`: Dataframe obtained executing the `data_preprocessing.R` file. It contains the molecular descriptors of 9213 chemical compounds after data preprocessing.
- `df_smiles.csv`: Dataframe obtained using *extract_smiles* function from `get_molecular_descriptors.R` file. It contains the smiles of 9213 chemical compounds.
