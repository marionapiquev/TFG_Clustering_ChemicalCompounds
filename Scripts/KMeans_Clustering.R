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
#LOAD THE DATAFRAMEs
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
}

########################################
#DATA PREPARETION (PCA)

#Prepare the data for clustering
DF <- Df_scaled #Choose the dataframe you want to use
cnames <- colnames(DF)
df <- DF[,cnames[!cnames %in% c("X","id")]]

#Perform PCA and choose the number of principle components
pca_DF <- pca(df)
pca_DF <- pca_DF[,0:2]
pca_DF['id'] <- DF$id
p <- plot_ly(pca_DF, x=pca_DF$Comp.1, y=pca_DF$Comp.2, 
             z=pca_DF$Comp.3, text = DF$id) %>% add_markers(size=1.5)
p

#We have found that there is an outlier, we remove this molecule
outlier <- 'DB11588'
pca_DF <- pca_DF %>%filter(id!=outlier)
scaled_DF <- Df_scaled %>%  filter(id!=outlier)
p <- plot_ly(pca_DF, x=pca_DF$Comp.1, y=pca_DF$Comp.2, 
             z=pca_DF$Comp.3) %>% add_markers(size=1.5)
p
########################################
#K-MEANS
#Prepare the data for clustering
DF <- pca_DF #Choose the dataframe you want to use
cnames <- colnames(DF)
df <- DF[,cnames[!cnames %in% c("X","id")]]

#Choose the number of clusters
set.seed(123)
bss <- c()
for (i in 1:50) {bss <- c(bss, kmeans(df,centers=i)$betweenss)}
plot(1:50, bss, type="l", xlab="Número de clústers",ylab="Suma de quadrats entre grupos")

wi <- kmeans(df,centers=1)$tot.withinss
for (i in 2:50) {wi[i] <- kmeans(df,centers=i)$tot.withinss}
plot(1:50, wi, type="l", xlab="Número de clústers",ylab="Suma de quadrats entre grupos")


#We calculate the different indexes for the range of k
metrics <- data.frame(k = integer(), Dunn = double(), silhouette = double())
dist_matrix <- daisy(df, metric = "euclidean")
d <- dist(dist_matrix)
for (i in 4:9){
  km.res <- kmeans(df, i)
  info <- cluster.stats(d,km.res$cluster)
  metrics[nrow(metrics) + 1,] <- data.frame(k = i, info$dunn2, silhouette = info$avg.silwidth)
}
#write.csv(metrics, "metrics.csv")

#Perform the clustering with the choosed number of k
km.res <- kmeans(df, 5)
km.res

#Get Dunn Index and Silhouette Index
info <- cluster.stats(d,km.res$cluster)
cat('Dunn Index:', info$dunn2)
cat('Shilouette Index:', info$avg.silwidth)

#Add clustering results to a dataframe using PCA dataframe
cnames <- colnames(scaled_DF)
df <- scaled_DF[,cnames[!cnames %in% c("X","id")]]
pca_DF <- pca(df)
pca_DF <- pca_DF[,0:3]

df_km <- pca_DF
df_km['cluster'] <- km.res['cluster']

#Visualize the clustering
fig <- plot_ly(df_km, x = ~Comp.1, y = ~Comp.2, z = ~Comp.3, color = ~cluster)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'C1'),yaxis = list(title = 'C2'),zaxis = list(title = 'C3')))
fig

#Save the final Dataframe with results
df_km <- scaled_DF
df_km['cluster'] <- km.res['cluster']
write.csv(df_km, "df_km.csv")

########################################
#CLUSTER ANALYSIS

#Look for the descriptors that influence more on deciding the label
cnames<- colnames(df_km)
cor_matrix <- cor(df_km[,cnames[!cnames %in% c("X","id","cluster")]],df_km$cluster, method = "pearson")
colnames(cor_matrix) <- 'cluster'
cor_df <- data.frame(cor_matrix) %>% filter(cluster > 0.4 | cluster < -0.4)
ggcorrplot(cor_df, colors = c("lavender", "steelblue", "navyblue"),lab = TRUE)

#Visulaize the boxplots
par(mfrow = c(3, 3))
for(i in rownames(cor_df)){
  boxplot(as.double(unlist(df_km[i])) ~ df_km$cluster, xlab='Cluster', ylab =i,col=c("lavender", "slategray2", "skyblue3", "steelblue", "navyblue"))
}


