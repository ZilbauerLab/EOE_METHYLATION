# Check for optimal number of clusters for each dataset
# 24Nov20 (fp215)
# 
# Folders required:
#     EOE_METHYLATION/CLUSTER_PLOTS/OPTIMAL_CLUSTERS
# 
# Files required:
#     EOE_METHYLATION/OUTPUT_DATASETS/EoE_NormalisedFiltered_GRSet.RData (containing normalised, QCed dataset produced by 01_EoE_RGsetNormalisationQC.v2.6.R)


#######################
# WORKING ENVIRONMENT #
#######################

setwd("~/EOE_METHYLATION")

library(plyr)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(ggplot2)
library(ClassDiscovery)
library(factoextra)
library(clustertend)
library(NbClust)
library(rafalib)
library(dendextend)
source("EoE_tools.R")

load("OUTPUT_DATASETS/EoE_NormalisedFiltered_GRSet.RData")

readme.epic.grset

samples$T0.eos <- samples$eos
samples$T0.eos[samples$dg.tpt=="EoE.T1"] <- NA
samples$T0.eos[samples$dg.tpt=="EoE.T2"] <- NA
samples$T1.eos <- samples$eos
samples$T1.eos[samples$dg.tpt=="EoE.T0"] <- NA
samples$T1.eos[samples$dg.tpt=="EoE.T2"] <- NA

# Convert beta values to M values for analysis:

convert <- function(x){log2(x/(1-x))}

EPIC.M <- structure(sapply(beta, convert), dim=dim(beta))
rownames(EPIC.M) <- rownames(beta)
colnames(EPIC.M) <- colnames(beta)

# Exclude outlier (SZ1.T0), T2 sample (SZ4.T2) & biologicical duplicates:

exclusions <- c("SZ1.T0","SZ4.T2","SZ1.T1.b","SZ19.b","SZ4.T0.b","SZ5.T1.b","SZ9.b")

all.samples <- subset(samples, !(sample.id %in% exclusions))
T0.samples <- subset(all.samples, dg.tpt!="EoE.T1")
T1.samples <- subset(all.samples, dg.tpt!="EoE.T0")
eoe.samples <- subset(all.samples, diagnosis=="EoE")

T0.beta <- beta[,T0.samples$array.id]
T1.beta <- beta[,T1.samples$array.id]
eoe.beta <- beta[,eoe.samples$array.id]
all.beta <- beta[,all.samples$array.id]

T0.M <- EPIC.M[,T0.samples$array.id]
T1.M <- EPIC.M[,T1.samples$array.id]
eoe.M <- EPIC.M[,eoe.samples$array.id]
all.M <- EPIC.M[,all.samples$array.id]

# Save analysis ready datasets for future use:

readme.subsets <- data.frame(dataset=c("beta","EPIC.M","samples","all.beta","all.M","all.samples","T0.beta","T0.M","T0.samples","T1.beta","T1.M","T1.samples","eoe.beta","eoe.M","eoe.samples"),description=c("Funnorm normalised, filtered EoE EPIC methylation array beta values (n=31)","Funnorm normalised, filtered EPIC methylation array M values (n=31)","Updated phenotype information for EPIC.beta and EPIC.M","EPIC.beta with SZ1.T0, SZ4.T2 & duplicates removed (n=24)","EPIC.M with SZ1.T0, SZ4.T2 & duplicates removed (n=24)","Updated phenotype information for all.beta and all.M","EPIC.beta subset for EoE pre-treatment T0 and control samples, with SZ1.T0 & duplicates removed (n=19)","EPIC.M subset for EoE pre-treatment T0 and control samples, with SZ1.T0 & duplicates removed (n=19)","Updated phenotype information for T0.EPIC and T0.M","EPIC.beta subset for EoE post-treatment T1 and control samples, duplicates removed (n=18)","EPIC.M subset for EoE post-treatment T1 and control samples, duplicates removed (n=18)","Updated phenotype information for T1.EPIC and T1.M","EPIC.beta subset for EoE T0 and T1 samples only, with SZ1.T0 & duplicates removed (n=11)","EPIC.M subset for EoE T0 and T1 samples only, with SZ1.T0 & duplicates removed (n=11)","Updated phenotype information for eoe.EPIC and eoe.M"))
save(beta,EPIC.M,samples,all.beta,all.M,all.samples,T0.beta,T0.M,T0.samples,T1.beta,T1.M,T1.samples,eoe.beta,eoe.M,eoe.samples,readme.subsets,file="OUTPUT_DATASETS/EoE_NormalisedFiltered_Subsets.RData")


#########################################################
# DETERMINE OPTIMAL NUMBER OF CLUSTERS FOR EACH DATASET #
#########################################################

# Plot using kmeans, PAM & hcut methods for comparison

require(gridExtra)

T0.filt <- T0.beta[!rowSums(!is.finite(T0.beta)),]
colnames(T0.filt) <- T0.samples$dups
T1.filt <- T1.beta[!rowSums(!is.finite(T1.beta)),]
colnames(T1.filt) <- T1.samples$dups
eoe.filt <- eoe.beta[!rowSums(!is.finite(eoe.beta)),]
colnames(eoe.filt) <- eoe.samples$dups

# Controls vs EoE.T0

pdf("CLUSTER_PLOTS/OPTIMAL_CLUSTERS/EoE-T0_vs_Control_OptimalClusters.pdf", width=15, height=5)
p1 <- fviz_nbclust(t(T0.filt), kmeans, method="silhouette")
p2 <- fviz_nbclust(t(T0.filt), pam, method="silhouette")
p3 <- fviz_nbclust(t(T0.filt), hcut, method="silhouette")
grid.arrange(p1+labs(subtitle="kmeans"),p2+labs(subtitle="PAM"),p3+labs(subtitle="hcut"),ncol=3)
dev.off()

# Controls vs EoE.T1

pdf("CLUSTER_PLOTS/OPTIMAL_CLUSTERS/EoE-T1_vs_Control_OptimalClusters.pdf", width=15, height=5)
p1 <- fviz_nbclust(t(T1.filt), kmeans, method="silhouette")
p2 <- fviz_nbclust(t(T1.filt), pam, method="silhouette")
p3 <- fviz_nbclust(t(T1.filt), hcut, method="silhouette")
grid.arrange(p1+labs(subtitle="kmeans"),p2+labs(subtitle="PAM"),p3+labs(subtitle="hcut"),ncol=3)
dev.off()

EoE.T0 vs EoE.T1

pdf("CLUSTER_PLOTS/OPTIMAL_CLUSTERS/EoE-T0_vs_EoE-T1_OptimalClusters.pdf", width=15, height=5)
p1 <- fviz_nbclust(t(eoe.filt), kmeans, method="silhouette")
p2 <- fviz_nbclust(t(eoe.filt), pam, method="silhouette")
p3 <- fviz_nbclust(t(eoe.filt), hcut, method="silhouette")
grid.arrange(p1+labs(subtitle="kmeans"),p2+labs(subtitle="PAM"),p3+labs(subtitle="hcut"),ncol=3)
dev.off()


################
# SESSION INFO #
################

sessionInfo()
