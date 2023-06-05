########################################
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
#LOADING R PACKAGES
library(rcdk)
library(ChemmineR) 
library(ChemmineOB) 
library(rpubchem)
library(chemblr)
library(fingerprint)
library(Rcpi)
library(cluster)
library("factoextra")
library(fpc)
library(ggcorrplot)
library(plotly)
library(ggplot2)
library(dbscan)
library(matrixStats)
library(dplyr)
library(RxnSim)
########################################
#LOAD THE DATAFRAME
#Set the directory 
setwd("/Users/MarionaP/Desktop/TFG2")
Df_scaled <- read.csv("df_molecular_descriptors_scaled.csv")
df_original <- read.csv("df_molecular_descriptors.csv")
smiles <- read.csv("df_smiles.csv")
########################################
#FUNCTIONS

#Function that returns the PCA 
pca <- function(df){
  pca.DF <- princomp(df) 
  print(summary(pca.DF))
  pca_DF <- as.data.frame(pca.DF$scores)
  return(pca_DF)
}

#Function that computes Tanimoto coeficient using a similarity matrix.
mols_similarity_matrix <- function(df, cnames, nclusters){
  for (k in 1:nclusters){
    cluster <- filter(df, cluster == k)
    c_id <- cluster[cnames[cnames %in% c("id")]]
    c_smiles <- filter(smiles, id %in% as.list(c_id)$id)
    m <- ms.compute.sim.matrix(c_smiles$smile, format = 'SMILES', fp.type='pubchem', sim.method='tanimoto')
    #print(m)
    m[lower.tri(m, diag=TRUE)] <- 0
    similarity <- sum(m)/sum(rowCounts(m > 0))
    cat('\nTanimoto coef mean cluster',k,':',similarity)
  }
}

########################################
#DATA PREPARATION (PCA)

#Prepare the data for clustering
DF <- Df_scaled #Choose the dataframe you want to use
cnames <- colnames(DF)
df <- DF[,cnames[!cnames %in% c("X","id")]]

#Perform PCA and choose the number of principle components
pca_DF <- pca(df)
pca_DF <- pca_DF[,0:2] #Choose number of components
pca_DF['id'] <- DF$id
p <- plot_ly(pca_DF, x=pca_DF$Comp.1, y=pca_DF$Comp.2, 
             z=pca_DF$Comp.3, text = DF$id) %>% add_markers(size=1.5) 
p <- p %>% layout(scene = list(xaxis = list(title = 'PC1'),yaxis = list(title = 'PC2'),zaxis = list(title = 'PC3')))
p

#We have found an outlier (DB11588), we remove this molecule
outlier <- 'DB11588'
pca_DF <- pca_DF %>%  filter(id!=outlier)
scaled_DF <- Df_scaled %>%  filter(id!=outlier)
df_original <- df_original %>%  filter(id!=outlier)
p <- plot_ly(pca_DF, x=pca_DF$Comp.1, y=pca_DF$Comp.2, 
             z=pca_DF$Comp.3) %>% add_markers(size=1.5)
p
########################################
#HIERARCHICAL CLUSTERING
#Prepare the data for HC
DF <- pca_DF #Choose the dataframe you want to use
cnames <- colnames(DF)
df <- DF[,cnames[!cnames %in% c("X","id")]]
dist_matrix <- daisy(df, metric = "euclidean")

#Perform the clustering
hc <- hclust(dist_matrix, method = 'ward.D')
plot(hc,main = paste('Dendrogram'), xlab = 'Molecules', ylab = 'Euclidean distances')

#Choose the optimal number of clusters
metrics <- data.frame(k = integer(), Dunn = double(), silhouette = double())
for (i in 3:12){
  grp <- cutree(hc, k = i)
  info <- cluster.stats(dist_matrix,grp)
  metrics[nrow(metrics) + 1,] <- data.frame(k = i, info$dunn2, silhouette = info$avg.silwidth)
}
#write.csv(metrics, "metrics.csv")

#Cut tree with the choosed number of k
cut <- 3
grp <- cutree(hc, k = cut)
plot(hc,main = paste('Dendrogram'),xlab = 'Molecules', ylab = 'Euclidean distances') 
rect.hclust(hc, k = cut, border = 2:5)

#Add clustering results to a dataframe using PCA dataframe
cnames <- colnames(scaled_DF)
df <- scaled_DF[,cnames[!cnames %in% c("X","id")]]
pca_DF <- pca(df)
pca_DF <- pca_DF[,0:3]
df_hc <- pca_DF
df_hc['cluster'] <- grp

#Visualize the clustering
fig <- plot_ly(df_hc, x = ~Comp.1, y = ~Comp.2, z = ~Comp.3, color = ~as.character(cluster), legendgroup=~cluster,  colors = c("slategray2","steelblue", "skyblue3"))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'C1'),yaxis = list(title = 'C2'),zaxis = list(title = 'C3')))
fig


#Add clustering results to a dataframe using scaled dataframe
df_hc <- scaled_DF
df_hc['cluster'] <- grp
#write.csv(df_hc, "df_hc.csv")

#Get Dunn Index and Silhouette Index
info <- cluster.stats(dist_matrix,df_hc$cluster)
cat('Dunn Index:', info$dunn2)
cat('Shilouette Index:', info$avg.silwidth)

#Get Cluster Tanimoto Similarity 
mols_similarity_matrix(df_hc, colnames(df_hc), cut)

########################################
#CLUSTER ANALYSIS
df_hc <- df_original
df_hc['cluster'] <- grp

#Look for the descriptors that influence more on deciding the label
cnames<- colnames(df_hc)
cor_matrix <- cor(df_hc[,cnames[!cnames %in% c("X","id","cluster")]],df_hc$cluster, method = "pearson")
colnames(cor_matrix) <- 'cluster'
cor_df <- data.frame(cor_matrix) %>% filter(cluster > 0.4 | cluster < -0.4)
ggcorrplot(cor_df, colors = c("lavender", "steelblue", "navyblue"),lab = TRUE)

#Visulaize the boxplots
par(mfrow = c(3, 3))
for(i in rownames(cor_df)){
  boxplot(as.double(unlist(df_hc[i])) ~ df_hc$cluster, xlab='Cluster', ylab =i,col=c("lavender", "slategray2", "skyblue3", "steelblue", "navyblue"))
}
