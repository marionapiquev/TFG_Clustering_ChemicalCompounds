########################################
#LOADING R PACKAGES
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
library(rcdk)
library(ChemmineR) 
library(ChemmineOB) 
library(rpubchem)
library(chemblr)
library(Rcpi)
library(cluster)
library("factoextra")
library(fpc)
library(ggcorrplot)
library(plotly)
library(ggplot2)
library(dbscan)
library(matrixStats)
library(RxnSim)
library(hrbrthemes)
library(RColorBrewer)
library(dplyr)
########################################
#LOAD THE DATAFRAME
#Set the directory 
setwd("/Users/MarionaP/Desktop/TFG2")
Df_scaled <- read.csv("df_molecular_descriptors_scaled.csv")
df_original <- read.csv("df_molecular_descriptors.csv")
smiles <- read.csv("df_smiles.csv")
df_filter <- read.csv("df_molecular_descriptors_scaled_filtered.csv")
########################################
#FUNCTIONS
#Function to visualize molecules for each cluster
visualize_cluster_mols <- function(df_dbscan, cnames, nclusters){
  for (i in 1:nclusters){
    cluster1 <- filter(df_dbscan, cluster == i)
    c1_id <- cluster1[cnames[cnames %in% c("id")]]
    c1_smiles <- filter(smiles, id %in% as.list(c1_id)$id)
    c1_sdf <- c()
    for (i in as.list(c1_smiles['smile'])$smile){
      c1_sdf <- c(c1_sdf, smiles2sdf(i))
    }
    par(mfrow=c(4,4))
    for (i in 1:length(c1_sdf)){
      openBabelPlot(c1_sdf[[i]])
    }
  }
}

#Function to get the similarity comparison of molecules on the same cluster.
#It computes the Tanimoto coeficient and the Overlap Coeficient.
mols_similarity <- function(df_dbscan, cnames, nclusters){
  mol_similarity <- data.frame(mol1 =  character(), mol2 =  character(), Tanimoto_Coef = double(), Overlap_Coef = double(), cluster = integer())
  for (k in 1:nclusters){  
    cluster1 <- filter(df_dbscan, cluster == k)
    c1_id <- cluster1[cnames[cnames %in% c("id")]]
    c1_smiles <- filter(smiles, id %in% as.list(c1_id)$id)
    c_filter <- c1_smiles
    for (i in c1_smiles['smile']$smile){
      c_filter <- filter(c_filter, smile != i)
      for (j in c_filter['smile']$smile){
        sim <- calcDrugMCSSim(i, j, type = 'smile')
        mol_similarity[nrow(mol_similarity) + 1,] <- data.frame(mol1 = i, mol2 = j, Tanimoto_Coef = sim[[2]], Overlap_Coef = sim[[3]], cluster = k)
      }
    }
  }
  for (i in 1:nclusters){
    cluster <- filter(mol_similarity, cluster == i)
    cat('\nTanimoto coef mean cluster',i,':',sum(cluster['Tanimoto_Coef'])/nrow(cluster))
    cat('\nOverlap coef mean cluster',i,':',sum(cluster['Overlap_Coef'])/nrow(cluster))
  }
  return(mol_similarity)
}

#Function that computes Tanimoto coeficient using a similarity matrix.
mols_similarity_matrix <- function(df_dbscan, cnames, nclusters){
  for (k in 1:nclusters){
    cluster <- filter(df_dbscan, cluster == k)
    c_id <- cluster[cnames[cnames %in% c("id")]]
    c_smiles <- filter(smiles, id %in% as.list(c_id)$id)
    m <- ms.compute.sim.matrix(c_smiles$smile, format = 'SMILES', fp.type='pubchem', sim.method='tanimoto')
    #print(m)
    m[lower.tri(m, diag=TRUE)] <- 0
    similarity <- sum(m)/sum(rowCounts(m > 0))
    cat('\nTanimoto coef mean cluster',k,':',similarity)
  }
}

#Function that returns the PCA 
pca <- function(df){
  pca.DF <- princomp(df) 
  print(summary(pca.DF))
  pca_DF <- as.data.frame(pca.DF$scores)
  return(pca_DF)
}

#Function that calculates the performance of the DBSCAN algorithm for different parameters values based on the Tanimoto Similarity
parameters_dbscan <- function(df, eps, minPts){
  mol_similarity <- data.frame(eps = double(), minPts = integer(), nclusters = integer(), non_class = integer(), cluster_max_size = integer(), Similarity = double())
  cnames <- colnames(df_original)
  for (i in eps){
    for (j in minPts){ 
      print(i);print(j)
      total <- 0
      h_dbscan_result <- dbscan(df, eps = i,minPts = j)
      df_hdbscan <- df_original
      df_hdbscan['cluster'] <- h_dbscan_result$cluster
      no_class <- nrow(filter(df_hdbscan, cluster == 0))
      max_size <- filter(df_hdbscan, cluster > 0) %>% dplyr::count(cluster)
      c_max_size <- max(max_size['n'])
      print(c_max_size)
      nclusters <- length(unique(df_hdbscan$cluster))-1
      if (nclusters > 0){
        for (k in 1:nclusters){
          cluster <- filter(df_hdbscan, cluster == k)
          if (nrow(cluster) < 3000 & nrow(cluster) > 1){
            c_id <- cluster[cnames[cnames %in% c("id")]]
            c_smiles <- filter(smiles, id %in% as.list(c_id)$id)
            m <- ms.compute.sim.matrix(c_smiles$smile, format = 'SMILES', fp.type='pubchem', sim.method='tanimoto')
            m[lower.tri(m, diag=TRUE)] <- 0
            similarity <- sum(m)/sum(rowCounts(m > 0))
            if(similarity == 'NaN'){total <- total + 0 } #if similarity is 0 the result is NaN
            else{total <- total + similarity}
          }
          else{
            total <- total + 0
          }
        }
        mol_similarity[nrow(mol_similarity) + 1,] <- data.frame(eps = i, minPts = j, nclusters = (length(unique(h_dbscan_result$cluster))-1), non_class = no_class, cluster_max_size = c_max_size, Similarity = total/(length(unique(h_dbscan_result$cluster))-1))
      }
      else{
        mol_similarity[nrow(mol_similarity) + 1,] <- data.frame(eps = i, minPts = j, nclusters = (length(unique(h_dbscan_result$cluster))-1), non_class = no_class, cluster_max_size = c_max_size, Similarity = 0)
      }
    }
  }
  return(mol_similarity)
}

########################################
#DATA PREPARATION (PCA)
#Prepare the data for clustering
DF <- Df_scaled #Choose the dataframe you want to use
cnames <- colnames(DF)
df <- DF[,cnames[!cnames %in% c("X","id")]]

pca_DF <- pca(df)
pca_DF <- pca_DF[,0:2] #Choose number of principle components
p <- plot_ly(pca_DF, x=pca_DF$Comp.1, y=pca_DF$Comp.2, 
             z=pca_DF$Comp.3) %>% add_markers(size=1.5)
print(p)
########################################
#DBSCAN
#Prepare the data for clustering
DF <-  Df_scaled #Choose the dataframe you want to use
cnames <- colnames(DF)
df <- DF[,cnames[!cnames %in% c("X","id")]]

#Calculate the performance of the algorithm with different eps and minPts based on the similarity between molecules on the same cluster
eps <- seq(2, 10, by=1)
minPts <- seq(2, 12, by=1)
param_similarity <- parameters_dbscan(df, eps, minPts)
write.csv(param_similarity, "param.csv")

mycolors <- colorRampPalette(brewer.pal(9, "Blues"))(15)
plot_dbscan <- ggplot(param_similarity, aes(eps, Similarity, colour = factor(minPts))) +
  geom_point()+ scale_fill_manual(mycolors)+geom_line()+theme_ipsum() + ggtitle("Mean Tanimoto Similarity") 
plot_dbscan
plot_dbscan_clusters <- ggplot(param_similarity, aes(eps, nclusters, colour = factor(minPts))) +
  geom_point()+ scale_fill_manual(mycolors)+geom_line()+theme_ipsum() + ggtitle("Number of clusters") 
plot_dbscan_clusters

#Perform the clustering
dbscan_result <- dbscan(df,eps = 2, minPts = 5)

#Add clustering results to a dataframe
df_dbscan <- df_original
df_dbscan['cluster'] <- dbscan_result$cluster

#Visualize the clustering plot
df_clusters <- filter(df_dbscan, cluster > 0)
ggplot(data = df_clusters, mapping = aes(x = MLogP,y = OB_MW, color = cluster))+geom_point() 

#fig <- plot_ly(df_clusters, x = ~MLogP, y = ~OB_MW, z = ~OB_logP, color = ~cluster)
#fig <- fig %>% add_markers()
#fig <- fig %>% layout(scene = list(xaxis = list(title = 'MLogP'),yaxis = list(title = 'OB_MW'),zaxis = list(title = 'OB_logP')))
#fig

#Visulaize the molecules of each cluster
nclusters <- length(unique(dbscan_result$cluster))-1
visualize_cluster_mols(df_dbscan, cnames,nclusters)

#Compute the similarity between molecules
#Tanimoto and Overlap coeficients (If the size of the clusters is big use the other function mols_similarity_matrix())
df_mols_similarity <- mols_similarity(df_dbscan, cnames, nclusters)

#Tanimoto coeficient
mols_similarity_matrix(df_dbscan, cnames, nclusters)


