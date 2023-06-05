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
#IMPORT MOLECULES IN R 
#Set the directory
setwd("/Users/MarionaP/Desktop/TFG2")

#Load molecules using a Database
mols <- load.molecules("DrugBank_3Dstructures.sdf")

#Load dataframes for results and visualitzation
Df_scaled <- read.csv("df_molecular_descriptors_scaled.csv")
df_original <- read.csv("df_molecular_descriptors.csv")
smiles <- read.csv("df_smiles.csv")
########################################
#FUNCTIONS

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

outlier <- 'DB11588'
########################################
#GET FINGERPRINTS
#We extract the fingerprints of the compunds
fps <- lapply(mols, get.fingerprint, type='pubchem') #Choose type

########################################
#PREPROCESSING
#Compute the intermolecular distances by the Tanimoto index
fp.sim <- fp.sim.matrix(fps, method='tanimoto')
fp.dist <- 1 - fp.sim

#MDS
fit <- cmdscale(dist(fp.dist), k=2)
colnames(fit) <- c("Dim.1", "Dim.2")
df_mds <- as.data.frame(fit)
#write.csv(df_mds, "mds.csv")
p <- plot_ly(mds, x=df_mds$Dim.1, y=df_mds$Dim.2) %>% add_markers(size=1.5)
########################################
#K-MEANS
#Choose the number of clusters
set.seed(123)
bss <- c()
for (i in 1:50) {bss <- c(bss, kmeans(fp.dist,centers=i)$betweenss)}
plot(1:50, bss, type="l", xlab="Número de clústers",ylab="Suma de quadrats entre grupos")

wi <- kmeans(fp.dist,centers=1)$tot.withinss
for (i in 2:50) {wi[i] <- kmeans(fp.dist,centers=i)$tot.withinss}
plot(1:50, wi, type="l", xlab="Número de clústers",ylab="Suma de quadrats entre grupos")

#We calculate the different indexes for the range of k
metrics <- data.frame(k = integer(), Dunn = double(), silhouette = double())
d <- dist(fp.dist)
for (i in 4:12){
  km.res <- kmeans(d, i)
  info <- cluster.stats(d,km.res$cluster)
  metrics[nrow(metrics) + 1,] <- data.frame(k = i, info$dunn2, silhouette = info$avg.silwidth)
}
#write.csv(metrics, "metrics.csv")

#Perform the clustering with the choosed number of k
centroids <- 4
d <- dist(fp.dist) 
km.res <- kmeans(d, centroids)
km.res

#Save the final Dataframe with results
df_km <- Df_scaled
df_km['cluster'] <- km.res['cluster']
write.csv(df_km, "df_km_finger.csv")

#Add clustering results to a dataframe using PCA dataframe
DF <- Df_scaled #Choose the dataframe you want to use
cnames <- colnames(DF)
df <- DF[,cnames[!cnames %in% c("X","id")]]
pca_DF <- pca(df)[,0:3]
pca_DF['id'] <- Df_scaled$id
df_km <- pca_DF 
df_km['cluster'] <- km.res['cluster']
df_km <- df_km%>%filter(id!=outlier)

#Visualize the clustering
fig <- plot_ly(df_km, x = ~Comp.1, y = ~Comp.2, z = ~Comp.3, color = ~as.character(cluster), legendgroup=~cluster,  colors = c("slategray2", "dodgerblue4"))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'C1'),yaxis = list(title = 'C2'),zaxis = list(title = 'C3')))
fig

#Get Dunn Index and Silhouette Index
dist_matrix <- daisy(df, metric = "euclidean")
d <- dist(dist_matrix)
info <- cluster.stats(d,km.res$cluster)
cat('Dunn Index:', info$dunn2)
cat('Shilouette Index:', info$avg.silwidth)

#Get Cluster Tanimoto Similarity 
mols_similarity_matrix(df_km, colnames(df_km), centroids)

p <- plot_ly(df, x=df$Dim.1, y=df$Dim.2, color = as.character(km.res$cluster),  colors = c("skyblue2", "steelblue","slategray2","navyblue")) %>% add_markers(size=1.5)
p <- plot_ly(df_mds, x=df_mds$Dim.1, y=df_mds$Dim.2, z=df_mds$Dim.3,color = as.character(df_km$cluster)) %>% add_markers(size=1.5)
########################################
#K-MEANS + MDS
mds <- read.csv("mds_extended.csv")
cnames <- colnames(mds)
df <- mds[,cnames[!cnames %in% c("X")]]

#Choose the number of clusters
set.seed(123)
bss <- c()
for (i in 1:50) {bss <- c(bss, kmeans(df,centers=i)$betweenss)}
plot(1:50, bss, type="l", xlab="Número de clústers",ylab="Suma de quadrats entre grupos")

wi <- kmeans(df,centers=1)$tot.withinss
for (i in 2:50) {wi[i] <- kmeans(df,centers=i)$tot.withinss}
plot(1:50, wi, type="l", xlab="Número de clústers",ylab="Suma de quadrats entre grupos")

#We calculate the different indexes for the range of k
metrics <- data.frame(k = integer(),  withinss= double(), silhouette = double())
dist_matrix <- daisy(df, metric = "euclidean")
d <- dist(dist_matrix)
for (i in 3:12){
  km.res <- kmeans(df, i)
  info <- cluster.stats(d,km.res$cluster)
  metrics[nrow(metrics) + 1,] <- data.frame(k = i, km.res$tot.withinss, silhouette = info$avg.silwidth)
}
#write.csv(metrics, "metrics.csv")

#Perform the clustering with the choosed number of k
centroids <- 4
km.res <- kmeans(df, centroids)
km.res

#Save the final Dataframe with results
df_km <- Df_scaled
df_km['cluster'] <- km.res['cluster']
#write.csv(df_km, "df_km_finger.csv")

#Add clustering results to a dataframe using PCA dataframe
DF <- Df_scaled #Choose the dataframe you want to use
cnames <- colnames(DF)
df2 <- DF[,cnames[!cnames %in% c("X","id")]]
pca_DF <- pca(df2)[,0:3]
pca_DF['id'] <- Df_scaled$id
df_km <- pca_DF 
df_km['cluster'] <- km.res['cluster']
#df_km <- df_km%>%filter(id!=outlier)

#Visualize the clustering
fig <- plot_ly(df_km, x = ~Comp.1, y = ~Comp.2, z = ~Comp.3, color = ~as.character(cluster), legendgroup=~cluster,  colors = c("slategray2", "dodgerblue4"))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'C1'),yaxis = list(title = 'C2'),zaxis = list(title = 'C3')))
fig

#Get Dunn Index and Silhouette Index
info <- cluster.stats(d,df_km$cluster)
cat('Dunn Index:', info$dunn2)
cat('Shilouette Index:', info$avg.silwidth)

#Get Cluster Tanimoto Similarity 
mds <- read.csv("mds.csv")
cnames <- colnames(mds)
df <- mds[,cnames[!cnames %in% c("X")]]
mols_similarity_matrix(df_km, colnames(df_km), centroids)

p <- plot_ly(df, x=df$Dim.1, y=df$Dim.2, color = as.character(km.res$cluster),  colors = c("skyblue2", "steelblue","slategray2","navyblue")) %>% add_markers(size=1.5)
########################################
#Hierarchical Clustering
#Prepare the data for HC
d <- dist(fp.dist) 
distance <- as.dist(d) 

#Perform the clustering
hc <- hclust(distance, method = 'ward.D')
plot(hc,main = paste('Dendrogram'), xlab = 'Molecules', ylab = 'Euclidean distances')

#Choose the optimal number of clusters between the range 4 a 10
metrics <- data.frame(k = integer(), Dunn = double(), silhouette = double())
for (i in 4:20){
  grp <- cutree(hc, k = i)
  info <- cluster.stats(distance,grp)
  metrics[nrow(metrics) + 1,] <- data.frame(k = i, info$dunn2, silhouette = info$avg.silwidth)
}
#write.csv(metrics, "metrics.csv")

#Cut tree with the choosed number of k
cut <- 3
grp <- cutree(hc, cut)
plot(hc,main = paste('Dendrogram'),xlab = 'Molecules', ylab = 'Euclidean distances') 
rect.hclust(hc, k = cut, border = 2:5)

#Save the final Dataframe with results
df_hc <- Df_scaled
df_hc['cluster'] <- grp
#write.csv(df_km, "df_hc_finger.csv")

#Add clustering results to a dataframe using PCA dataframe
DF <- Df_scaled #Choose the dataframe you want to use
cnames <- colnames(DF)
df <- DF[,cnames[!cnames %in% c("X","id")]]
pca_DF <- pca(df)[,0:3]
pca_DF['id'] <- Df_scaled$id
df_hc <- pca_DF 
df_hc['cluster'] <- grp
df_hc <- df_hc%>%filter(id!=outlier)

#Visualize the clustering
fig <- plot_ly(df_hc, x = ~Comp.1, y = ~Comp.2, z = ~Comp.3, color = ~as.character(cluster), legendgroup=~cluster,  colors = c("slategray2", "dodgerblue4"))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'C1'),yaxis = list(title = 'C2'),zaxis = list(title = 'C3')))
fig

#Get Dunn Index and Silhouette Index
info <- cluster.stats(d,grp)
cat('Dunn Index:', info$dunn2)
cat('Shilouette Index:', info$avg.silwidth)

#Get Cluster Tanimoto Similarity 
mols_similarity_matrix(df_hc, colnames(df_hc), cut)

########################################
#Hierarchical Clustering + MDS
#Prepare the data for HC
mds <- read.csv("mds.csv")
cnames <- colnames(mds)
df <- mds[,cnames[!cnames %in% c("X")]]
dist_matrix <- daisy(df, metric = "euclidean")

#Perform the clustering
hc <- hclust(dist_matrix, method = 'average')
plot(hc,main = paste('Dendrogram'), xlab = 'Molecules', ylab = 'Euclidean distances')

#Choose the optimal number of clusters between the range 4 a 10
metrics <- data.frame(k = integer(), Dunn = double(), silhouette = double())
for (i in 3:15){
  grp <- cutree(hc, k = i)
  info <- cluster.stats(dist(dist_matrix),grp)
  metrics[nrow(metrics) + 1,] <- data.frame(k = i, info$dunn2, silhouette = info$avg.silwidth)
}
#write.csv(metrics, "metrics.csv")

#Cut tree with the choosed number of k
cut <- 4
grp <- cutree(hc, cut)
plot(hc,main = paste('Dendrogram'),xlab = 'Molecules', ylab = 'Euclidean distances') 
rect.hclust(hc, k = cut, border = 2:5)

#Save the final Dataframe with results
df_hc <- Df_scaled
df_hc['cluster'] <- grp
#write.csv(df_km, "df_hc_finger.csv")

#Add clustering results to a dataframe using PCA dataframe
DF <- Df_scaled #Choose the dataframe you want to use
cnames <- colnames(DF)
df2 <- DF[,cnames[!cnames %in% c("X","id")]]
pca_DF <- pca(df2)[,0:3]
pca_DF['id'] <- Df_scaled$id
df_hc <- pca_DF 
df_hc['cluster'] <- grp
#df_hc <- df_hc%>%filter(id!=outlier)
#Add clustering results to a dataframe using PCA dataframe

#df_km <- df_km%>%filter(id!=outlier)
#Visualize the clustering
fig <- plot_ly(df_hc, x = ~Comp.1, y = ~Comp.2, z = ~Comp.3, color = ~as.character(cluster), legendgroup=~cluster,  colors = c("slategray2", "dodgerblue4"))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'C1'),yaxis = list(title = 'C2'),zaxis = list(title = 'C3')))
fig

#Get Dunn Index and Silhouette Index
info <- cluster.stats(d,grp)
cat('Dunn Index:', info$dunn2)
cat('Shilouette Index:', info$avg.silwidth)

#Get Cluster Tanimoto Similarity 
mols_similarity_matrix(df_hc, colnames(df_hc), cut)

p <- plot_ly(df, x=df$Dim.1, y=df$Dim.2, color = as.character(grp),  colors = c("slategray2","steelblue", "navyblue")) %>% add_markers(size=1.5)

########################################
#CLUSTER ANALYSIS
#Get heatmap of the Tanimoto Similiarity inside each cluster
for (k in 1:nclusters){
  cluster <- filter(df_hc, cluster == k)
  cnames <- colnames(cluster)
  c_id <- cluster[cnames[cnames %in% c("id")]]
  c_smiles <- filter(smiles, id %in% as.list(c_id)$id)
  m <- ms.compute.sim.matrix(c_smiles$smile, format = 'SMILES', fp.type='pubchem', sim.method='tanimoto')
  pheatmap(as.matrix(m))
}


#Look for the descriptors that influence more on deciding the label
df_hc <- df_original
df_hc['nAtom'] <- extractDrugAtomCount(mols, silent = TRUE)
df_hc['nBond'] <-extractDrugBondCount(mols, silent = TRUE)
df_hc['cluster'] <- grp

cnames<- colnames(df_hc)
cor_matrix <- cor(df_hc[,cnames[!cnames %in% c("X","id","cluster")]],df_hc$cluster, method = "pearson")
colnames(cor_matrix) <- 'cluster'
cor_df <- data.frame(cor_matrix) %>% filter(cluster > 0.3 | cluster < -0.3)
ggcorrplot(cor_df, colors = c("lavender", "steelblue", "navyblue"),lab = TRUE)

par(mfrow = c(3, 3))
for(i in rownames(cor_df)){
  boxplot(as.double(unlist(df_hc[i])) ~ df_hc$cluster, xlab='Cluster', ylab =i,col=c("lavender", "slategray2", "skyblue3", "steelblue", "navyblue"))
}
