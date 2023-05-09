# RGset creation, filtering, normalisation and initial QC
# 18Nov20 (fp215)
#
# Folders required:
#     EOE_METHYLATION
#     EOE_METHYLATION/OUTPUT_DATASETS
#     EOE_METHYLATION/QC_PLOTS
#     EOE_METHYLATION/QC_PLOTS/01_RAW_DATA
#     EOE_METHYLATION/QC_PLOTS/02_NORMALISED_DATA
#     EOE_METHYLATION/QC_PLOTS/03_NORMALISED_FILTERED_DATA
#
# Files required:
#     EOE_METHYLATION/RAW_IDATS (containing the raw iDat files)
#     EOE_METHYLATION/EoE_phenotypes.txt (containing relevant phenotype & clinical info)


#######################
# WORKING ENVIRONMENT #
#######################

setwd("~/EOE_METHYLATION")

baseDir <- "~/EOE_METHYLATION/RAW_IDATS/"

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
library(quantro)
source("~/EoE_tools.R")
source("~/quantro2.R")


#################
# LOAD RAW DATA #
#################

# Info of interest: sample.id, case.no, diagnosis, timepoint, sex, age, eos (eosinophils per hpf), array.id, sentrix.id (batch)

samples <- read.table("EoE_phenotypes.txt", header=T, sep="\t", stringsAsFactors=F)
samples$timepoint[samples$diagnosis=="Control"] <- "T0"
samples$dg.tpt <- paste0(samlpes$diagnosis,".",samples$timepoint)
samples$dg.tpt[samples$dg.tpt=="Control.T0"] <- "Control"

RGset.raw <- read.metharray(file.path(baseDir, samples$array.id), force=T)
RGset.raw@annotation <- c(array="IlluminaHumanMethylationEPIC", annotation="ilm10b4.hg19")

pData(RGset.raw)$case.no <- samples$case.no
pData(RGset.raw)$array.id <- samples$array.id
pData(RGset.raw)$diagnosis <- samples$diagnosis
pData(RGset.raw)$timepoint <- samples$timepoint
pData(RGset.raw)$dg.tpt <- samples$dg.tpt
pData(RGset.raw)$sex <- samples$sex
pData(RGset.raw)$age <- samples$age
pData(RGset.raw)$eos <- samples$eos
pData(RGset.raw)$sentrix <- samples$sentrix.id

validObject(RGset.raw)


############
# PROBE QC #
############

# Calculate probes failing detection threshold using detectionP
# Remove samples with a high proportion of failed probes

pval <- detectionP(RGset.raw)
pval.mean <- colMeans(detectionP(RGset.raw))
samples$mean.detectionP <- pval.mean
pData(RGset.raw)$mean.detectionP <- pval.mean

pdf("QC_PLOTS/01_RAW_DATA/EoE_RawData_MeanDetectionP.pdf")
ggplot(samples)+geom_boxplot(aes(as.factor(sentrix), pval.mean, fill=as.factor(sentrix)), outlier.shape=NA)+geom_point(aes(as.factor(sentrix), pval.mean, group=array.id, fill=as.factor(sentrix)), shape=21, color="black", position=position_jitter(w=0.25))+theme_bw()+theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))+xlab("Sentrix ID")+ylab("Mean Detection P Value")+guides(fill=F)+ylim(0,0.008)
dev.off()

keep <- colMeans(pval) < 0.005
samples <- subset(samples, mean.detectionP < 0.005)
RGset <- RGset.raw[,keep]

identical(colnames(RGset),samples$array.id)

# Produce QC report for raw data filtered for poor quality samples:

qcReport(RGset, pdf="QC_PLOTS/01_RAW_DATA/EoE_RawData_QCReport.pdf")


###########
# QUANTRO #
###########

# Check whether quantile normalisation is needed and if there is a difference between sample groups
# Original quantro source code edited to exclude NAs, original code here: https://github.com/stephaniehicks/quantro/blob/master/R/AllClasses.R

beta.raw <- getBeta(RGset)
qtestPerm <- quantro(object=beta.raw[complete.cases(beta.raw),], groupFactor=samples$dg.tpt, B=1000)

qtestPerm

pdf("QC_PLOTS/01_RAW_DATA/EoE_QuantroPermutationTest_Diagnosis-Timepoint.pdf")
quantroPlot(qtestPerm)
dev.off()


##################
# NORMALISE DATA #
##################

# Normalise using preprocessFunnorm: https://www.rdocumentation.org/packages/minfi/versions/1.18.4/topics/preprocessFunnorm

GRset <- preprocessFunnorm(RGset, nPCs=2, bgCorr=T, dyeCorr=T)
beta.raw <- getBeta(RGset)
beta <- getBeta(GRset)

# Check beta distribution

set.seed(1)

beta.plot <- melt(beta[sample(nrow(beta),nrow(beta)/4)),])
beta.plot <- beta.plot[which(!(is.na(beta.plot$value))),]
colnames(beta.plot) <- c("CpG","ID","Beta")
beta.plot <- as.data.frame(merge(beta.plot, samples, by.x="ID", by.y="array.id"))

pdf("QC_PLOTS/02_NORMALISED_DATA/EoE_Normalised_BetaDistribution_DiagnosisTimepoint.pdf")
ggplot(beta.plot, aes(Beta, group=as.character(ID), color=as.character(dg.tpt)))+geom_density()+theme_bw()+colscale_eoe+xlab("DNAm Beta Value")
dev.off()

beta.filt <- beta[!rowSums(!is.finite(beta)),]

pdf("QC_PLOTS/02_NORMALISED_DATA/EoE_Normalised_BySampleBetaDistribution_DiagnosisTimepoint.pdf")
matboxplot(beta.filt, groupFactor=samples$dg.tpt, xaxt="n", main="Beta Value")
dev.off()

# PCA of first 10 PCs, prior to probe filtering and batch correction

pca <- prcomp(t(beta.filt))
pca <- as.data.frame(pca$x)
pca$array.id <- rownames(pca)
pca <- join(pca, samples, type="left", match="all")
write.table(pca, "QC_PLOTS/02_NORMALISED_DATA/EoE_Normalised_PCs.txt", row.names=F, sep="\t", quote=F)

pdf("QC_PLOTS/02_NORMALISED_DATA/EoE_Normalised_PCA_DiagnosisTimepoint.pdf")
splom(~pca[c("PC1","PC2","PC3","PC4","PC5")], groups=pca$dg.tpt, col=c("grey75","darkblue","darkorchid","darkred"), pch=rep(19,4), pscales=0, key=list(space="bottom", title="PCA: all normalised samples (beta values)", points=list(pch=rep(19,4), col=c("grey75","darkblue","darkorchid","darkred")), text=list(c("Control","EoE T0","EoE T1","EoE T2"))))
dev.off()

pdf("QC_PLOTS/02_NORMALISED_DATA/EoE_Normalised_PCA_Chip.pdf")
splom(~pca[c("PC1","PC2","PC3","PC4","PC5")], groups=pca$sentrix, col=c("#F8766D","#7CAE00","#00BFC4","#C77CFF"), pch=rep(19,4), pscales=0, key=list(space="bottom", title="PCA: all normalised samples (beta values)", points=list(pch=rep(19,4), col=c("#F8766D","#7CAE00","#00BFC4","#C77CFF")), text=list(c("203845510034","203845520020","203845520066","203845520080"))))
dev.off()


############
# SNP DATA #
############

# Extract genotyping probe data for reference

SNPSs <- getSnpBeta(RGset)
SNPs <- SNPs[complete.cases(SNPs),]
d <- dist(t(SNPs))
hc <- hclust(d, method="complete")
write.table(SNPs, "OUTPUT_DATASETS/EoE_GenotypingData.txt", row.names=F, sep="\t", quote=F)

pdf("QC_PLOTS/02_NORMALISED_DATA/EoE_RawData_ClusteringBySNP.pdf", width=30, height=15)
myplclust(hc, labels=samples$array.id, lab.col=as.fumeric(as.character(samples$case.no)), cex=1)
dev.off()


###################
# PROBE FILTERING #
###################

# Remove sex chromosomes, probes overlapping SNPs, cross-hybridising to other locations, with over 5% NAs or with a detectionP > 0.05 in 1% or more of samples

Fset <- dropLociWithSnps(GRset, snps=c("SBE","CpG"), maf=0)
beta <- getBeta(Fset)
beta <- rmSNPandCH(beta, dist=2, mafcut=0.05, rmcrosshyb=T, rmXY=T)

probes.all <- nrow(RGset)
probes.mapped <- nrow(GRset)
probes.nosnp <- nrow(Fset)
probes.nosnpch <- nrow(beta)

na.count <- sapply(1:nrow(beta), function(y) length(which(is.na(beta[y,]))))
na.count <- which(na.count < (ncol(beta)*0.05))
beta <- beta[na.count,]

probes.nona <- nrow(beta)

fails <- pval > 0.05
fail.p <- names(which(rowMeans(fails) > 0.01))
fail.id <- names(which(colMeans(fails) > 0.01))

beta <- beta[which(!(rownames(beta) %in% fail.p)),]
beta <- beta[,which(!(colnames(beta) %in% fail.id))]

probes.pval <- nrow(beta)

# Summarise probe filtering:

summary <- data.frame(Data=c("Raw RGset probes","Mapped probes","Probes overlapping SNPs removed","Cross-hyb & sex chr probes removed","Probes with over 5% NA removed","Probes with high detectionP in 1% samples removed"),Probe.Number=c(probes.all,probes.mapped,probes.nosnp,probes.nosnpch,probes.nona,probes.pval))
summary$Removed <- c(0,sapply(2:nrow(summary), function(x) summary$Probe.Number[x-1] - summary$Probe.Number[x]))
write.table(summary, "QC_PLOTS/03_NORMALISED_FILTERED_DATA/EoE_EPICData_ProbeFilteringSummary.txt", row.names=F, quote=F, sep="\t")
summary$Data <- factor(summary$Data, rev(summary$Data))

pdf("QC_PLOTS/03_NORMALISED_FILTERED_DATA/EoE_ProbeAttritionPlot.pdf", width=12)
ggplot(summary)+geom_bar(aes(Data,-Probe.Number), stat="identity", fill="grey70", color="black")+geom_bar(aes(Data,Removed), stat="identity", fill="darkred", color="black")+geom_text(aes(x=Data, y=-min(Probe.Number)/2, label=comma(Probe.Number)))+geom_text(aes(x=Data, y=max(Removed)/1.5, label=comma(Removed)))+geom_hline(yintercept=0)+coord_flip()+theme_bw()+ylab("")+xlab("")+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(colour="grey20",size=12),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_blank(),panel.background=element_blank())+scale_x_discrete(position="top")
dev.off()

# Check beta distribution

set.seed(1)

beta.plot <- melt(beta[sample(nrow(beta),(nrow(beta)/4)),])
beta.plot <- beta.plot[which(!(is.na(beta.plot$value))),]
colnames(beta.plot) <- c("CpG","ID","Beta")
beta.plot <- as.data.frame(merge(beta.plot, samples, by.x="ID", by.y="array.id"))

pdf("QC_PLOTS/03_NORMALISED_FILTERED_DATA/EoE_NormalisedFiltered_BetaDistribution_DiagnosisTimepoint.pdf")
ggplot(beta.plot, aes(Beta, group=as.character(ID), color=as.character(dg.tpt)))+geom_density()+theme_bw()+colscale_eoe+xlab("DNAm Beta Value")
dev.off()

pdf("QC_PLOTS/03_NORMALISED_FILTERED_DATA/EoE_NormalisedFiltered_BySampleBetaDistribution_DiagnosisTimepoint.pdf")
matboxplot(beta, groupFactor=samples$dg.tpt, xaxt="n", main="Beta Value")
dev.off()


####################################
# PCA OF NORMALISED, FILTERED DATA #
####################################

samples$dups <- ifelse(samples$diagnosis=="Control",samples$case.no,(paste0(samples$case.no,".",samples$timepoint)))
pca <- prcomp(t(beta))
pca <- as.data.frame(pca$x)
pca$array.id <- rownames(pca)
pca2 <- pca[1:10]
pca2$array.id <- pca$array.id
pca <- join(pca, samples, type="left", match="all")
write.table(pca, "QC_PLOTS/03_NORMALISED_FILTERED_DATA/EoE_NormalisedFiltered_PCs.txt", row.names=F, sep="\t", quote=F)

pdf("QC_PLOTS/03_NORMALISED_FILTERED_DATA/EoE_NormalisedFiltered_PCA_DiagnosisTimepoint.pdf")
splom(~pca[c("PC1","PC2","PC3","PC4","PC5")], groups=pca$dg.tpt, col=c("grey75","darkblue","darkorchid","darkred"), pch=rep(19,4), pscales=0, key=list(space="bottom", title="PCA: all normalised, filtered samples (beta values)", points=list(pch=rep(19,4), col=c("grey75","darkblue","darkorchid","darkred")), text=list(c("Control","EoE T0","EoE T1","EoE T2"))))
dev.off()

pdf("QC_PLOTS/03_NORMALISED_FILTERED_DATA/EoE_NormalisedFiltered_PCA_Chip.pdf")
splom(~pca[c("PC1","PC2","PC3","PC4","PC5")], groups=pca$sentrix, col=c("#F8766D","#7CAE00","#00BFC4","#C77CFF"), pch=rep(19,4), pscales=0, key=list(space="bottom", title="PCA: all normalised samples (beta values)", points=list(pch=rep(19,4), col=c("#F8766D","#7CAE00","#00BFC4","#C77CFF")), text=list(c("203845510034","203845520020","203845520066","203845520080"))))
dev.off()

pdf("QC_PLOTS/03_NORMALISED_FILTERED_DATA/EoE_NormalisedFiltered_PC1vPC2_DiagnosisTimepoint.pdf", width=6, height=4)
ggplot(pca, aes(PC1, PC2, fill=dg.tpt))+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+labs(colour="Diagnosis-Timepoint")+fillscale_eoe
dev.off()

pdf("QC_PLOTS/03_NORMALISED_FILTERED_DATA/EoE_NormalisedFiltered_PC1vPC2_Chip.pdf", width=6, height=4)
ggplot(pca, aes(PC1, PC2, fill=sentrix))+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual("Chip",values=c("#F8766D","#7CAE00","#00BFC4","#C77CFF"))
dev.off()

pdf("QC_PLOTS/03_NORMALISED_FILTERED_DATA/EoE_NormalisedFiltered_PC1vPC2_EosinophilsPerHPF.pdf", width=6, height=4)
ggplot(pca, aes(PC1, PC2, fill=eos.2))+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual("Eosinophils",values=c("red","blue"))
dev.off()

                              
###########################################
# CLUSTERING OF NORMALISED, FILTERED DATA #
###########################################

samples$col.dg.tpt <- as.factor(samples$dg.tpt)
samples$col.dg.tpt <- factor(samples$col.dg.tpt, levels=c("Control","EoE.T0","EoE.T1","EoE.T2"))
levels(samples$col.dg.tpt) <- c("grey75","darkblue","darkorchid","darkred")
samples$col.dg.tpt <- as.character(samples$col.dg.tpt)
                              
beta.filt <- beta[!rowSums(!is.finite(beta)),]
d <- dist(t(beta.filt))
hc <- hculst(d, method="complete")

# Plot clustering, identifying case.no and colour coding by combined diagnosis and timepoint:

pdf("QC_PLOTS/03_NORMALISED_FILTERED_DATA/EoE_NormalisedFiltered_Clustering_DiagnosisTimepoint.pdf", width=30, height=10)
myplclust(hc, labels=samples$case.no, lab.col=samples$col.dg.tpt, cex=2, main="Sample Clustering (complete), showing diagnosis & treatment timepoint")
legend("topright",legend=c("Control","EoE T0","EoE T1","EoE T2"),text.col=c("grey75","darkblue","darkorchid","darkred"),pch=rep(16,4),col=c("grey75","darkblue","darkorchid","darkred"))
dev.off()


#################
# SAVE DATASETS #
#################

# Raw data

readme.rgset <- data.frame(Dataset-c("RGset.raw","RGset","samples"),Description=c("RGChannelSet of raw EPIC array intensity data","RGset.raw filtered for samples with a high proportion of failed probes (mean.detectionP < 0.005)","Sample Info"))
save(RGset.raw,RGset,samples,readme.rgset,file="OUTPUT_DATASETS/EoE_RawData_RGSet.RData")

# Normalised, filtered data
                              
EPIC.GRset <- makeGenomicRatioSetFromMatrix(beta, array="IlluminaHumanMethylationEPIC", annotation="ilm10b4.hg19", pData=samples, mergeManifest=T, what="Beta")
             
EPIC.GRset

readme.epic.grset <- data.frame(Dataset=c("EPIC.GRset","beta","samples"),Description=c("Funnorm normalised, filtered EPIC GenomicRatioSet, poor quality samples and probes removed","Beta values for EPIC.GRset","Phenotype table for EPIC.GRset"))
save(EPIC.GRset,beta,samples,readme.epic.grset,file="OUTPUT_DATASETS/EoE_NormalisedFiltered_GRSet.RData")


################
# SESSION INFO #
################

sessionInfo()
