########################################
#INCREASE THE SIZE OF JAVA RUNTIME HEAP
options(java.parameters = "-Xmx")
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
########################################
#IMPORT MOLECULES IN R 

#Set the directory
setwd("/Users/MarionaP/Desktop/TFG2")

#Load molecules using a Database
#mols <- load.molecules("DrugBank_3Dstructures.sdf")
#sdfset <- read.SDFset("DrugBank_3Dstructures.sdf") # 9213 molecules

########################################
#EXTRACT PROPIETIES OF THE COMPOUNDS

#Extract molecular descriptors
extract_propietites_compounds <- function(name) {
  #LOAD MOLECULES 
  mols <- load.molecules(name)
  sdfset <- read.SDFset(name)
  #OPEN BABLE DESCRIPTORS
  smiles <- datablocktag(sdfset, tag="SMILES")
  ids <- datablocktag(sdfset, tag="DRUGBANK_ID")
  info <- extractDrugDescOB(smiles[1], type = 'smile')
  OB_structural <- data.frame(id=ids[1],OB_HBA1=info$HBA1, OB_HBA2=info$HBA2, OB_HBD=info$HBD,OB_nF=info$nF,OB_MW=info$MW)
  OB_physicochemical <- data.frame(id=ids[1],OB_MR=info$MR, OB_logP=info$logP, OB_TPSA=info$TPSA)
  for (i in 2:length(smiles)){
    info <- extractDrugDescOB(smiles[i], type = 'smile')
    OB_structural[nrow(OB_structural) + 1,] <- data.frame(id=ids[i],OB_HBA1=info$HBA1, OB_HBA2=info$HBA2, OB_HBD=info$HBD,OB_nF=info$nF,OB_MW=info$MW)
    OB_physicochemical[nrow(OB_physicochemical) + 1,] <- data.frame(id=ids[i],OB_MR=info$MR, OB_logP=info$logP, OB_TPSA=info$TPSA)
  }
  
  #STRUCTURAL DESCRIPTORS
  structural_descriptors <- data.frame(extractDrugCarbonTypes(mols, silent = TRUE),
                                       extractDrugAromaticAtomsCount(mols, silent = TRUE),
                                       extractDrugRotatableBondsCount(mols, silent = TRUE),
                                       extractDrugLargestChain(mols, silent = TRUE),
                                       extractDrugLargestPiSystem(mols, silent = TRUE),
                                       extractDrugLongestAliphaticChain(mols, silent = TRUE))
  
  structural_descriptors['id'] <- datablocktag(sdfset, tag="DRUGBANK_ID")
  #structural_descriptors['HybRatio'] <- extractDrugIPMolecularLearning(mols, silent = TRUE)
  
  #PHYSICOCHEMICAL DESCRIPTORS
  physicochemical_descriptors <- data.frame(extractDrugMannholdLogP(mols, silent = TRUE))
  extractDrugCPSA_df <- extractDrugCPSA(mols, silent = TRUE)
  physicochemical_descriptors['THSA'] <- extractDrugCPSA_df['THSA']
  physicochemical_descriptors['TPSA'] <- extractDrugCPSA_df['TPSA']
  physicochemical_descriptors['RPSA'] <- extractDrugCPSA_df['RPSA']
  physicochemical_descriptors['id'] <- datablocktag(sdfset, tag="DRUGBANK_ID")
  
  #GEOMETRICAL DESCRIPTORS
  geometrical_descriptors <- data.frame(extractDrugVABC(mols, silent = TRUE),
                                        extractDrugECI(mols, silent = TRUE),
                                        extractDrugFMF(mols, silent = TRUE),
                                        extractDrugVAdjMa(mols, silent = TRUE),
                                        extractDrugFragmentComplexity(mols, silent = TRUE),
                                        extractDrugZagrebIndex(mols, silent = TRUE))
  
  geometrical_descriptors['id'] <- datablocktag(sdfset, tag="DRUGBANK_ID")
  #geometrical_descriptors['MolIP'] <- extractDrugIPMolecularLearning(mols, silent = TRUE)
  
  #ELECTRONIC DESCRIPTORS
  electronic_descriptors <- extractDrugCPSA_df[,1:15]
  electronic_descriptors['RPCG'] <- extractDrugCPSA_df['RPCG']
  electronic_descriptors['RNCG'] <- extractDrugCPSA_df['RNCG']
  electronic_descriptors['id'] <- datablocktag(sdfset, tag="DRUGBANK_ID")
  
  #MERGE THE OPEN BABLE DATAFRAMES WITH THE STRUCTURAL AND PHYSIOCHEMICAL
  df_structural_descriptors <- merge(structural_descriptors,OB_structural,by='id')
  df_physicochemical_descriptors <- merge(physicochemical_descriptors,OB_physicochemical,by='id')
  
  #MERGE ALL THE DATAFRAMES
  df_molecular_descriptors1 <- merge(geometrical_descriptors,electronic_descriptors,by='id')
  df_molecular_descriptors2 <- merge(df_structural_descriptors,df_physicochemical_descriptors,by='id')
  df_molecular_descriptors <- merge(df_molecular_descriptors1,df_molecular_descriptors2,by='id')
  
  #IF YOU WANT TO SAVE EACH DF
  write.csv(df_molecular_descriptors, "df_molecular_descriptors.csv")
  #write.csv(df_structural_descriptors, "structural_descriptors.csv") 
  #write.csv(df_physicochemical_descriptors, "physicochemical_descriptors.csv") 
  #write.csv(geometrical_descriptors, "geometrical_descriptors.csv") 
  #write.csv(electronic_descriptors, "electronic_descriptors.csv") 
  return(df_molecular_descriptors)
}

df_MD <- extract_propietites_compounds("DrugBank_3Dstructures.sdf")

#Extract fingerprints
extract_fingerprints <- function(name,ftype) {
  mols <- load.molecules("DrugBank_3Dstructures.sdf")
  if (ftype == 'pubchem'){
    fp <- lapply(mols, get.fingerprint, type='pubchem')
  }
  else if  (ftype == 'exteded'){
    fp <- lapply(mols, get.fingerprint, type='extended')
  }
  return(fp)
}

fp <- extract_fingerprints("DrugBank_3Dstructures.sdf", 'pubchem')

#Extract smiles
extract_smiles <- function(name){
  #LOAD MOLECULES 
  sdfset <- read.SDFset(name)
  #EXTRACT SMILES
  smiles <- datablocktag(sdfset, tag="SMILES")
  ids <- datablocktag(sdfset, tag="DRUGBANK_ID")
  df_smiles <- data.frame(id=ids, smile=smiles)
  write.csv(df_smiles, "df_smiles.csv")
}

smiles <- extract_smiles("DrugBank_3Dstructures.sdf")

