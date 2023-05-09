# Subset data for all samples to top 25 biologically significant DMPs and examine clustering
# 25Nov20 (fp215)
#
# Folders required:
#     EOE_METHYLATION/CLUSTERING_BY_DMP
#
# Files required:
#     EOE_METHYLATION/DMA_RESULTS/EoE_DifferentialMethylationAnalysisResults_DMPs.RData


#######################
# WORKING ENVIRONMENT #
#######################

library(plyr)
library(ggplot2)
library(reshape)
library(ClassDiscovery)
library(factoextra)
library(clustertend)
library(NbClust)
library(rafalib)
library(dendextend)
library(gridExtra)

setwd("~/EOE_METHYLATION")
load("DMA_RESULTS/EoE_DifferentialMethylationAnalysisResults_DMPs.RData")


###################
# CREATE DATASETS #
###################

data <- dma.results[with(dma.results, order(fdr.T0,abs(db.T0))),]
data <- data[1:25,]
cpgs <- data$CpG

all.samples <- subset(samples, sample.id!="SZ1.T0" & sample.id!="SZ4.T2")
all.beta <- EPIC.beta[,all.samples$array.id]
all.beta <- subset(all.beta, rownames(all.beta) %in% cpgs)

all(colnames(all.beta)==all.samples$array.id)

colnames(all.beta) <- all.samples$sample.id


##############################
# OPTIMAL NUMBER OF CLUSTERS #
##############################

pdf("CLUSTERING_BY_DMP/EoE-T0_vs_Control_Top25-DMPs_in_AllSamples_OptimalClusters_n29.pdf", width=15, height=5)
p1 <- fviz_nbclust(t(all.beta), kmeans, method="silhouette")
p2 <- fviz_nbclust(t(all.beta), pam, method="silhouette")
p3 <- fviz_nbclust(t(all.beta), hcut, method="silhouette")
grid.arrange(p1+labs(subtitle="kmeans"),p2+labs(subtitle="PAM"),p3+labs(subtitle="hcut"),ncol=3)
dev.off()


###################
# PLOT DENDROGRAM #
###################

all.bar <- c("grey75","darkblue","darkorchid")
all.bar <- all.bar[as.numeric(as.factor(all.samples$dg.tpt))]
all.bar2 <- c("gold","darkorange","red","grey95")
all.bar2 <- all.bar2[as.numeric(as.factor(all.samples$eos))]

d <- distanceMatrix(all.beta, "pearson")
hc <- hclust(d, method="average")
dend <- as.dendrogram(hc)
dend <- color_branches(dend, k=2, col=c("darkgreen","darkred"))

pdf("CLUSTERING_BY_DMP/EoE-T0_vs_Control_Top25-DMPs_in_AllSamples_n29.pdf", width=12, height=5)
par(mar=c(7,7,4,2))
plot(dend)
colored_bars(colors=cbind(all.bar2,all.bar), dend=dend, rowLabels=c("Eosinophils","Diagnosis/Timepoint"))
legend("topright", legend=c("Control","EoE at diagnosis","EoE after treatment","Under 15 per hpf","15-20 per hpf","20-50 per hpf","Over 50 per hpf"), pch=15, bty="n", col=c("grey75","darkblue","darkorchid","grey95","gold","darkorange","red"))
dev.off()


################
# SESSION INFO #
################

sessionInfo()
