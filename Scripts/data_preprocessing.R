########################################
#INSTALLING R PACKAGES (ONLY RUN ONCE!!!)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("ChemmineR")
#BiocManager::install("Rcpi")
#BiocManager::install("ChemmineOB")

#library(devtools)
#install_github("CDK-R/rpubchem", dependencies=TRUE)
#devtools::install_github('rajarshi/chemblr/package')
########################################
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
library(dplyr)
########################################
#LOAD THE DATAFRAME
#Set the directory 
setwd("/Users/MarionaP/Desktop/TFG/DrugBank")
DF <- read.csv("df_molecular_descriptors.csv")
smiles <- read.csv("df_smiles.csv")
########################################
#DATA PREPROCESSING
#Check the dimension of the dataframe
dim(DF)

#Count number of nans for each column
na_count <-sapply(DF, function(y) sum(length(which(is.na(y))))/sum(nrow(DF)))
na_count <- data.frame(na_count)

#For each column if the na_count is greater than 0.3 drop the column
df <- DF
for (i in colnames(df)){
  if (na_count[i,] > 0.3){
    df <- df[, !(names(df) %in% c(i))]
    }
}
dim(df)

#Drop rows with NA
df <- na.omit(df)
dim(df)

#Drop the repeated molecules
df <- df %>% distinct(id, .keep_all = TRUE)
df <- df[! colnames(df) %in% c("X")]
write.csv(df, "df_molecular_descriptors.csv")

smiles <- smiles %>% distinct(id, .keep_all = TRUE)
smiles <- smiles[! colnames(smiles) %in% c("X")]
write.csv(smiles, "df_smiles.csv")

#Summary: Check if the dataframe is scaled 
summary(df)

#Scale the numeric variables
cnames <- colnames(df)
df_scaled <- as.data.frame(scale(df[,-1]))
df_scaled$id <- df$id
#df_scaled <- cbind(df[,colnames_no], scaled_desc)

#Visulaize the box plots
par(mfrow = c(4, 7))
for (i in colnames(df_scaled)){
  if(i!='id'){
    boxplot(df_scaled[i], xlab=i, col=c("slategray2", "skyblue3", "steelblue"))
  }}

#Check id the dataframe is scaled
summary(df_scaled)

#Feature Selection
#Visualize correlation matrix
cor_matrix <- cor(df_scaled[,cnames[!cnames %in% c("X","id")]], method = "pearson")
ggcorrplot(cor_matrix, colors = c("lavender", "steelblue", "navyblue"),lab = FALSE)

#Drop the C2SP1 and C1SP1 variables since they don't correlate with the other variables (all correlations lower than 10%)
df_scaled <- df_scaled[, !(names(df_scaled) %in% c('C2SP1','C1SP1'))]
cnames <- colnames(df_scaled)
cor_matrix <- cor(df_scaled[,cnames[!cnames %in% c("X","id")]], method = "pearson")
ggcorrplot(cor_matrix, colors = c("lavender", "steelblue", "navyblue"),lab = FALSE)

#Download scaled dataframe
write.csv(df_scaled, "df_molecular_descriptors_scaled.csv")
