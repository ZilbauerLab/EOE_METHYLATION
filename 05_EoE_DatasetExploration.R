# Data exploration
# 19Nov20 (fp215)
# 
# Folders required:
#       EOE_METHYLATION/CLUSTER_PLOTS
#       EOE_METHYLATION/PCA_PLOTS
#       EOE_METHYLATION/HEAT_SCREE_PLOTS
#
# Files required:
#       EOE_METHYLATION/OUTPUT_DATASETS/EoE_NormalisedFiltered_Subsets.RData (containing normalised, QCed, subset datasets produced by 04_EoE_OptimalClusteringCheck.v1.1.R)


#######################
# WORKING ENVIRONMENT #
#######################

setwd("~/EOE_METHYLATION")

library(plyr)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(ggplot2)
library(lattice)
library(reshape)
library(scales)
library(ClassDiscovery)
library(factoextra)
library(clustertend)
library(NbClust)
library(rafalib)
library(dendextend)
source("EoE_tools.R")

load("OUTPUT_DATASETS/EoE_NormalisedFiltered_Subsets.RData")

readme.subsets


#################
# CLUSTER PLOTS #
#################

# NB, eos column: under.15 = < 15 eos/hpf, over.15 = 15 - 20 eos/hpf, over.20 = 20 - 50 eos/hpf, over.50 = > 50 eos/hpf
# Controls vs EoE.T0

all(colnames(T0.beta)==T0.samples$array.id)

T0.bar <- c("grey75","darkblue")
T0.bar <- T0.bar[as.numeric(as.factor(T0.samples$dg.tpt))]
T0.bar2 <- c("darkorange","red","grey95")
T0.bar2 <- T0.bar2[as.numeric(as.factor(T0.samples$T0.eos))]

beta.filt <- T0.beta[!rowSums(!is.finite(T0.beta)),]
colnames(beta.filt) <- T0.samples$case.no
d <- dist(t(beta.filt), "pearson")
hc <- hclust(d, method="average")
dend <- as.dendrogram(hc)
dend <- color_branches(dend, k=2, col=c("darkgreen","darkred"))

pdf("CLUSTER_PLOTS/EoE-T0_vs_Control_Clustering.pdf", width=12, height=6)
par(mar=c(7,7,4,2))
plot(dend)
colored_bars(colors=cbind(T0.bar2,T0.bar), dend=dend, rowLabels=c("Eosinophils","Diagnosis/Timepoint"))
legend("topright", legend=c("Control","EoE at diagnosis","Under 15 per hpf","20-50 per hpf","Over 50 per hpf"), pch=15, bty="n", col=c("grey75","darkblue","grey95","darkorange","red"))
dev.off()

# Controls vs EoE.T1

all(colnames(T1.beta)==T1.samples$array.id)

T1.bar <- c("grey75","darkorchid")
T1.bar <- T1.bar[as.numeric(as.factor(T1.samples$dg.tpt))]
T1.bar2 <- c("darkorange","grey95")
T1.bar2 <- T1.bar2[as.numeric(as.factor(T1.samples$T1.eos))]

beta.filt <- T1.beta[!rowSums(!is.finite(T1.beta)),]
colnames(beta.filt) <- T1.samples$case.no
d <- dist(t(beta.filt), "pearson")
hc <- hclust(d, method="average")
dend <- as.dendrogram(hc)
dend <- color_branches(dend, k=2, col=c("darkgreen","darkred"))

pdf("CLUSTER_PLOTS/EoE-T1_vs_Control_Clustering.pdf", width=12, height=6)
par(mar=c(7,7,4,2))
plot(dend)
colored_bars(colors=cbind(T1.bar2,T1.bar), dend=dend, rowLabels=c("Eosinophils","Diagnosis/Timepoint"))
legend("topright", legend=c("Control","EoE after treatment","Under 15 per hpf","20-50 per hpf"), pch=15, bty="n", col=c("grey75","darkorchid","grey95","darkorange"))
dev.off()

# EoE.T0 vs EoE.T1

all(colnames(eos.beta)==eos.samples$array.id)

eos.bar <- c("darkblue","darkorchid")
eos.bar <- eos.bar[as.numeric(as.factor(eos.samples$dg.tpt))]
eos.bar2 <- c("darkorange","red","grey95")
eos.bar2 <- eos.bar2[as.numeric(as.factor(eos.samples$eos))]

beta.filt <- eos.beta[!rowSums(!is.finite(eos.beta)),]
colnames(beta.filt) <- eos.samples$case.no
d <- dist(t(beta.filt), "pearson")
hc <- hclust(d, method="average")
dend <- as.dendrogram(hc)
dend <- color_branches(dend, k=2, col=c("darkgreen","darkred"))

pdf("CLUSTER_PLOTS/EoE-T0_vs_EoE-T1_Clustering.pdf", width=12, height=6)
par(mar=c(7,7,4,2))
plot(dend)
colored_bars(colors=cbind(eos.bar2,eos.bar), dend=dend, rowLabels=c("Eosinophils","Diagnosis/Timepoint"))
legend("topright", legend=c("EoE at diagnosis","EoE after treatment","Under 15 per hpf","20-50 per hpf","Over 50 per hpf"), pch=15, bty="n", col=c("darkblue","darkorchid","grey95","darkorange","red"))
dev.off()


#############
# PCA PLOTS #
#############

require(gridExtra)

# Controls vs EoE.T0

pca <- prcomp(t(T0.beta))
pca <- as.data.frame(pca$x)
pca$array.id <- rownames(pca)
pca2 <- pca[1:10]
pca2$array.id <- pca$array.id
pca <- join(pca, T0.samples, type="left", match="all")

pdf("PCA_PLOTS/EoE-T0_vs_Control_PC1vPC2.pdf", width=15, height=5)
p1 <- ggplot(pca, aes(PC1, PC2, fill=dg.tpt))
p2 <- ggplot(pca, aes(PC1, PC2, fill=T0.eos))
p3 <- ggplot(pca, aes(PC1, PC2, fill=sex))
grid.arrange((p1+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+labs(colour="Diagnosis-Timepoint")+fillscale_eoe),(p2+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual("T0 Eosinophils",values=c("darkorange","red","grey75"),labels=c("20-50 per hpc","Over 50 per hpc","Under 15 per hpc"))),(p3+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual("Gender",values=c("red","blue"),labels=c("Female","Male"))),ncol=3)
dev.off()

# Controls vs EoE.T1

pca <- prcomp(t(T1.beta))
pca <- as.data.frame(pca$x)
pca$array.id <- rownames(pca)
pca2 <- pca[1:10]
pca2$array.id <- pca$array.id
pca <- join(pca, T1.samples, type="left", match="all")

pdf("PCA_PLOTS/EoE-T1_vs_Control_PC1vPC2.pdf", width=15, height=5)
p1 <- ggplot(pca, aes(PC1, PC2, fill=dg.tpt))
p2 <- ggplot(pca, aes(PC1, PC2, fill=T1.eos))
p3 <- ggplot(pca, aes(PC1, PC2, fill=sex))
grid.arrange((p1+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+labs(colour="Diagnosis-Timepoint")+fillscale_eoe),(p2+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual("T1 Eosinophils",values=c("darkorange","grey75"),labels=c("20-50 per hpc","Under 15 per hpc"))),(p3+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual("Gender",values=c("red","blue"),labels=c("Female","Male"))),ncol=3)
dev.off()

# EoE.T0 vs EoE.T1

pca <- prcomp(t(eoe.beta))
pca <- as.data.frame(pca$x)
pca$array.id <- rownames(pca)
pca2 <- pca[1:10]
pca2$array.id <- pca$array.id
pca <- join(pca, eoe.samples, type="left", match="all")

pdf("PCA_PLOTS/EoE-T0_vs_EoE-T1_PC1vPC2.pdf", width=5, height=5)
ggplot(pca, aes(PC1, PC2, fill=dg.tpt))+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+labs(colour="Diagnosis-Timepoint")+fillscale_eoe)
dev.off()

# All samples (excluding outliers & duplicates)

pca <- prcomp(t(all.beta))
pca <- as.data.frame(pca$x)
pca$array.id <- rownames(pca)
pca2 <- pca[1:10]
pca2$array.id <- pca$array.id
pca <- join(pca, meta, type="left", match="all")
write.table(pca, "PCA_PLOTS/EoE_All-NoOutliersNoDups_PCs_n24.txt", row.names=F, sep="\t", quote=F)
pca$timepoint[pca$diagnosis=="Control"] <- "T0"

pdf("PCA_PLOTS/EoE_All-NoOutliersNoDups_PC1vPC2_Diagnosis-EosinophilsProgression_n24.pdf", width=9, height=7)
ggplot(pca, aes(PC1, PC2))+geom_point(aes(shape=as.factor(dg.tpt), color=as.factor(eos)), size=3)+geom_line(aes(PC1,PC2,group=case.no),color="black")+scale_shape_manual(name="Diagnosis & timepoint",values=c(19,15,17),labels=c("Control","EoE T0","EoE T1"))+scale_color_manual(name="Eosinophils",values=c("gold","darkorange","red","grey75"),labels=c("15-20 per hpf","20-50 per hpf","Over 50 per hpf","Under 15 per hpf"))+theme_bw()
dev.off()

pdf("PCA_PLOTS/EoE_All-NoOutliersNoDups_PC1-PC5-lattice_DiagnosisTimepoint_n24.pdf")
splom(~pca[c("PC1","PC2","PC3","PC4","PC5")], groups=pca$dg.tpt, col=c("grey75","darkblue","darkorchid"), pch=rep(19,3), pscales=0, key=list(space="bottom", title="PCA: all normalised, filtered samples (beta values)", points=list(pch=rep(19,3), col=c("grey75","darkblue","darkorchid")), text=list(c("Control","EoE T0","EoE T1"))))
dev.off()


####################
# HEAT SCREE PLOTS #
####################

PCs_to_view <- 10

# Controls vs EoE.T0

pca <- prcomp(t(T0.beta))
Loadings <- as.data.frame(pca$x)
vars <- pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- T0.samples[c("diagnosis","sex","sentrix","T0.eos")]
meta_continuous <- T0.samples[c("mean.detectionP","age")]
colnames(meta_categorical) <- c("Diagnosis","Sex","Sentrix ID","Eosinophils")
colnames(meta_continuous) <- c("Mean detectionP","Age")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("HEAT_SCREE_PLOTS/EoE-T0_vs_Control_HeatScreePlot.pdf", width=10, height=6)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()

# Controls vs EoE.T1

pca <- prcomp(t(T1.beta))
Loadings <- as.data.frame(pca$x)
vars <- pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- T1.samples[c("diagnosis","sex","sentrix","T1.eos")]
meta_continuous <- T1.samples[c("mean.detectionP","age")]
colnames(meta_categorical) <- c("Diagnosis","Sex","Sentrix ID","Eosinophils")
colnames(meta_continuous) <- c("Mean detectionP","Age")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("HEAT_SCREE_PLOTS/EoE-T1_vs_Control_HeatScreePlot.pdf", width=10, height=6)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()


################
# SESSION INFO #
################

sessionInfo()
