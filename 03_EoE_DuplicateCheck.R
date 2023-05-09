# Check correlation between biological duplicate samples
# 19Nov20 (fp215)
#
# Folders required:
#       EOE_METHYLATION/DUPLICATES
#
# Files required:
#       EOE_METHYLATION/OUTPUT_DATASETS/EoE_NormalisedFiltered_GRSet.RData (containing normalised, QCed dataset produced by 01_EoE_RGsetNormalisationQC.v2.6.R)
#       EOE_METHYLATION/03_NORMALISED_FILTERED_DATA/EoE_NormalisedFiltered_PCs.txt (sample info with first 10 PCs for each sample produced by 01_EoE_RGsetNormalisationQC.v2.6.R)


#######################
# WORKING ENVIRONMENT #
#######################

setwd("~/EOE_METHYLATION")

library(minfi)
library(DMRcate)
library(plyr)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(ggplot2)
library(lattice)
library(RColorBrewer)
library(rafalib)
library(reshape)
library(scales)
library(LSD)
library(ggpubr)
require(gridExtra)
source("EoE_tools.R")

load("OUTPUT_DATASETS/EoE_NormalisedFiltered_GRSet.RData")
pca <- read.table("03_NORMALISED_FILTERED_DATA/EoE_NormalisedFiltered_PCs.txt", header=T, sep="\t", stringsAsFactors=F)

readme.epic.grset


###############################
# EXTRACT DUPLICATES DATASETS #
###############################

dup.list <- c("SZ1.T1","SZ4.T0","SZ5.T1","SZ9.T0","SZ19.T0")
pca2 <- subset(pca, dups %in% dup.list)
dups <- subset(samples, dups %in% dup.list)

dups[c("case.no","timepoint","diagnosis","mean.detectionP","array.id","sample.id")]

beta <- as.data.frame(beta)
colnames(beta) <- samples$sample.id
beta <- beta[,dups$sample.id]
beta2 <- beta[sample(nrow(beta),(nrow(beta)/25)),]

readme.dups <- data.frame(Dataset=c("beta","dups"),Description=c("Funnorm normalised, filtered EPIC beta values for biological duplicates","Phenotype table for beta"))
save(beta,dups,readme.dups,file="OUTPUT_DATASETS/EoE_NormalisedFiltered_Duplicates.RData")

# Compare phenotypic info for biological duplicates:

dups[c("case.no","timepoint","diagnosis","mean.detectionP","array.id","sample.id")]


#########################
# DUPLICATE CORRELATION #
#########################

# Calculate correlation (Kendall):

cor.SZ1 <- cor(beta2$SZ1.T1.a, beta2$SZ1.T1.b, method="kendall")
cor.SZ4 <- cor(beta2$SZ4.T0.a, beta2$SZ4.T0.b, method="kendall")
cor.SZ5 <- cor(beta2$SZ5.T1.a, beta2$SZ5.T1.b, method="kendall")
cor.SZ9 <- cor(beta2$SZ9.a, beta2$SZ9.b, method="kendall")
cor.SZ19 <- cor(beta2$SZ19.a, beta2$SZ19.b, method="kendall")

# Calculate residual standard deviations:

lm.SZ1 <- sd(residuals(lm(beta2$SZ1.T1.a~beta2$SZ1.T1.b)))
lm.SZ4 <- sd(residuals(lm(beta2$SZ4.T0.a~beta2$SZ4.T0.b)))
lm.SZ5 <- sd(residuals(lm(beta2$SZ5.T1.a~beta2$SZ5.T1.b)))
lm.SZ9 <- sd(residuals(lm(beta2$SZ9.a~beta2$SZ9.b)))
lm.SZ19 <- sd(residuals(lm(beta2$SZ19.a~beta2$SZ19.b)))

# Create results table:

cor.results <- unique(dups[c("case.no","timepoint","diagnosis")])
cor.results$sample.1 <- c("SZ1.T1.a","SZ19.a","SZ4.T0.a","SZ5.T1.a","SZ9.a")
cor.results$sample.2 <- c("SZ1.T1.b","SZ19.b","SZ4.T0.b","SZ5.T1.b","SZ9.b")
cor.results$R <- c(cor.SZ1, cor.SZ19, cor.SZ4, cor.SZ5, cor.SZ9)
cor.results$R2 <- with(cor.results, R * R)
cor.results$SD.Residuals <- c(lm.SZ1, lm.SZ19, lm.SZ4, lm.SZ5, lm.SZ9)
write.table(cor.results, "DUPLICATES/EoE_BiologicalDuplicateCorrelation.txt", row.names=F, sep="\t", quote=F)

cor.results

# Plot results:

sz1 <- round((subset(cor.results, case.no=="SZ1")$R2),digits=2)
sz4 <- round((subset(cor.results, case.no=="SZ4")$R2),digits=2)
sz5 <- round((subset(cor.results, case.no=="SZ5")$R2),digits=2)
sz9 <- round((subset(cor.results, case.no=="SZ9")$R2),digits=2)
sz19 <- round((subset(cor.results, case.no=="SZ19")$R2),digits=2)

pdf("DUPLICATES/EoE_BiologicalDuplicateCorrelation.pdf", width=15, height=10)
p1 <- ggscatter(beta2, x="SZ1.T1.a", y="SZ1.T1.b", add="reg.line", conf.int=T, cor.coef=F, cor.method="kendall", xlab="Sample 1", ylab="Sample 2", title=paste0("SZ1.T1: R2=",sz1,", SD residuals= ",round(lm.SZ1,digits=2)))
p2 <- ggscatter(beta2, x="SZ4.T0.a", y="SZ4.T0.b", add="reg.line", conf.int=T, cor.coef=F, cor.method="kendall", xlab="Sample 1", ylab="Sample 2", title=paste0("SZ4.T0: R2=",sz4,", SD residuals= ",round(lm.SZ4,digits=2)))
p3 <- ggscatter(beta2, x="SZ5.T1.a", y="SZ5.T1.b", add="reg.line", conf.int=T, cor.coef=F, cor.method="kendall", xlab="Sample 1", ylab="Sample 2", title=paste0("SZ5.T1: R2=",sz5,", SD residuals= ",round(lm.SZ5,digits=2)))
p4 <- ggscatter(beta2, x="SZ9.a", y="SZ9.b", add="reg.line", conf.int=T, cor.coef=F, cor.method="kendall", xlab="Sample 1", ylab="Sample 2", title=paste0("SZ9: R2=",sz9,", SD residuals= ",round(lm.SZ9,digits=2)))
p5 <- ggscatter(beta2, x="SZ19.a", y="SZ19.b", add="reg.line", conf.int=T, cor.coef=F, cor.method="kendall", xlab="Sample 1", ylab="Sample 2", title=paste0("SZ19: R2=",sz19,", SD residuals= ",round(lm.SZ19,digits=2)))
grid.arrange(p1,p2,p3,p4,p5,ncol=3)
dev.off()


############################
# PCA OF DUPLICATE SAMPLES #
############################

pdf("DUPLICATES/EoE_BiologicalDuplicates_PC1vPC2.pdf", width=5, height=4)
ggplot(pca2, aes(PC1, PC2, fill=dups))+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_discrete(name="Replicates")
dev.off()


################
# SESSION INFO #
################

sessionInfo()
