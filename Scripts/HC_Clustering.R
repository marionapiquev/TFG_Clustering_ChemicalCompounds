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
pca <- function(df){
  pca.DF <- princomp(df) 
  print(summary(pca.DF))
  #Show the percentage of variances explained by each principal component.
  fviz_eig(pca.DF)
  #Graph of individuals (Individuals with a similar profile are grouped together)
  fviz_pca_ind(pca.DF,
               col.ind = "cos2", 
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE     # Avoid text overlapping
  )
  pca_DF <- as.data.frame(pca.DF$scores)
}

#Prepare the data for clustering
DF <- Df_scaled #Choose the dataframe you want to use
cnames <- colnames(DF)
df <- DF[,cnames[!cnames %in% c("X","id")]]

pca_DF <- pca(df)
pca_DF <- pca_DF[,0:2]
pca_DF['id'] <- DF$id
p <- plot_ly(pca_DF, x=pca_DF$Comp.1, y=pca_DF$Comp.2, 
             z=pca_DF$Comp.3, text = DF$id) %>% add_markers(size=1.5)
p

#We have found an outlier (DB11588), we remove this molecule
outlier <- 'DB11588'
pca_DF <- pca_DF %>%  filter(id!=outlier)
scaled_DF <- Df_scaled %>%  filter(id!=outlier)
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

#Choose the optimal number of clusters between the range 4 a 10
metrics <- data.frame(k = integer(), Dunn = double(), silhouette = double())
for (i in 3:12){
  grp <- cutree(hc, k = i)
  info <- cluster.stats(dist_matrix,grp)
  metrics[nrow(metrics) + 1,] <- data.frame(k = i, info$dunn2, silhouette = info$avg.silwidth)
}
write.csv(metrics, "metrics.csv")

#Cut tree with the choosed number of k
grp <- cutree(hc, k = 3)
plot(hc,main = paste('Dendrogram'),xlab = 'Molecules', ylab = 'Euclidean distances') 
rect.hclust(hc, k = 3, border = 2:5)

#Add clustering results to a dataframe using PCA dataframe
cnames <- colnames(scaled_DF)
df <- scaled_DF[,cnames[!cnames %in% c("X","id")]]
pca_DF <- pca(df)
pca_DF <- pca_DF[,0:3]

df_hc <- pca_DF
df_hc['cluster'] <- grp

#Visualize the clustering
fig <- plot_ly(df_hc, x = ~Comp.1, y = ~Comp.2, z = ~Comp.3, color = ~cluster)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'C1'),yaxis = list(title = 'C2'),zaxis = list(title = 'C3')))
fig

df_hc <- scaled_DF
df_hc['cluster'] <- grp
write.csv(df_hc, "df_hc.csv")

########################################
#CLUSTER ANALYSIS

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

