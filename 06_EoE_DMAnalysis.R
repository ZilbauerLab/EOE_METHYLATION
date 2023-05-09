# Differential methylation analysis
# 24Nov20 (fp215)
#
# Folders required:
#     EOE_METHYLATION/DMA_RESULTS
#
# Files required:
#     EOE_METHYLATION/OUTPUT_DATASETS/EoE_NormalisedFiltered_Subsets.RData (containing normalised, QCed, subset datasets produced by 04_EoE_OptimalClusteringCheck.v1.1.R)
#     Simplified_MethylationEPIC_v1-0_B4.txt (containing genomic info for Illumina EPIC CpGs - CpG, chromosome, GRCh37 position, Gene symbol etc)


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
library(rafalib)
library(reshape)
library(scales)
library(limma)
source("EoE_tools.R")
require(gridExtra)
set.seed(123)

info <- read.table("Simplified_MethylationEPIC_v1-0_B4.txt", header=T, sep="\t", stringsAsFactors=F)
load("OUTPUT_DATASETS/EoE_NormalisedFiltered_Subsets.RData")

readme.subsets


################
# DMA FUNCTION #
################

# Controls vs EoE

DMA <- function(mval,beta,meta){
  print(paste("Testing in ", ncol(beta), " individuals", sep=""))

  # Construct model and fit linear model & mixed model
  
  meta$dg.tpt <- as.factor(meta$dg.tpt)
  meta$dg.tpt <- factor(meta$dg.tpt, levels=c("Control","EoE.T0","EoE.T1"))
  meta$sex <- as.factor(meta$sex)
  meta$age <- as.numeric(meta$age)

  mod <- model.matrix(~0+dg.tpt+sex+age, data=meta)
  fit <- lmFit(mval, mod)
  corrfit <- duplicateCorrelation(mval, mod, block=meta$case.no)
  fit2 <- lmFit(mval, mod, block=meta$case.no, correlation=corrfit$consensus)
  
 # Construct the contrast matrices

  contrastMatrix <- makeContrasts(
    T0 = dg.tptControl - dg.tptEoE.T0,
    T1 = dg.tptControl - dg.tptEoE.T1,
    eoe = dg.tptEoE.T0 - dg.tptEoE.T1,
    levels = mod)
  
 # Compute contrasts for linear model fit and calculate differential expression using ebayes
  
  contrastFit <- contrasts.fit(fit, contrastMatrix)
  contrastFitEb <- eBayes(contrastFit)
  T0.contrast <- topTable(contrastFitEb, coef="T0", number=Inf)
  T1.contrast <- topTable(contrastFitEb, coef="T1", number=Inf)

  contrastFit2 <- contrasts.fit(fit2, contrastMatrix)
  contrastFitEb2 <- eBayes(contrastFit2)
  eoe.contrast <- topTable(contrastFitEb2, coef="eoe", number=Inf)

  # Calculate covariate adjusted beta values

  avebeta.lm <- apply(beta, 1, function(x){
    lm(x~sex+age, data=meta)})
  residuals <- t(sapply(avebeta.lm, function(x)residuals(summary(x))))
  colnames(residuals) <- colnames(beta)
  adj.residuals <- residuals+matrix(apply(beta, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))

  # Calculate delta beta

  delta_beta <- lapply(1:nrow(adj.residuals), function(x) {
    group_means <- tapply(adj.residuals[x,], meta$dg.tpt, mean)
    c(group_means[2]-group_means[1],group_means[3]-group_means[1],group_means[3]-group_means[2])
  })
  
  limmares <- as.data.frame(do.call(rbind, delta_beta))
  colnames(limmares) <- c("db.T0","db.T1","db.eoe")
  limmares$CpG <- rownames(adj.residuals)

  # Construct results table
  
  results.T0 <- data.frame(Pvalue.T0=T0.contrast$P.Value, logFC.T0=T0.contrast$logFC, CpG=rownames(T0.contrast))
  results.T1 <- data.frame(Pvalue.T1=T1.contrast$P.Value, logFC.T1=T1.contrast$logFC, CpG=rownames(T1.contrast))
  results.eoe <- data.frame(Pvalue.eoe=eoe.contrast$P.Value, logFC.eoe=eoe.contrast$logFC, CpG=rownames(eoe.contrast))
  
  res <- merge(results.T0, results.T1, by="CpG")
  res <- merge(res, results.eoe, by="CpG")
  limmares <- merge(limmares, res, by="CpG")

  limmares$fdr.T0 <- p.adjust(limmares$Pvalue.T0, method = "BH", n = nrow(limmares))
  limmares$fdr.T1 <- p.adjust(limmares$Pvalue.T1, method = "BH", n = nrow(limmares))
  limmares$fdr.eoe <- p.adjust(limmares$Pvalue.eoe, method = "BH", n = nrow(limmares))

  limmares <- join(limmares, info, type="left", match="all")
  
  limmares}


###########################################
# DIFFERENTIAL METHYLATION ANALYSIS: DMPs #
###########################################

dma.results <- DMA(all.M,all.beta,all.samples)

# Save results

write.table(dma.results, "DMA_RESULTS/EoE_DifferentialMethylationAnalysisResults_DMPs.txt", row.names=F, quote=F, sep="\t")

readme.dma <- data.frame(dataset=c("dma.results"),description=c("Limma DMA results for EoE, controls vs EoE.T0, controls vs EoE.T1 & EoE.T0 vs EoE.T1, excluding duplicates (n=24)"))
save(dma.results,readme.dma,file="DMA_RESULTS/EoE_DifferentialMethylationAnalysisResults_DMPs.RData")


############################
# CALCULATE NUMBER OF HITS #
############################

# Calculate the number of biologically relevant hits with the following criteria:
#       - delta beta (dB) >= 0.05
#       - fdr adjusted p-value <= 0.01
# NB, no results seen for EoE.T0 vs EoE.T1

hit_count <- function(dB, fdr, data){
  sig.fdr <- as.data.frame(data)[which(data$fdr.T0 <= fdr),]
  sig.dB <- sig.fdr[which(abs(sig.fdr$db.T0) >= dB),]
  print(paste("Control vs EoE.T0: ",nrow(sig.dB)))

  sig.fdr <- as.data.frame(data)[which(data$fdr.T1 <= fdr),]
  sig.dB <- sig.fdr[which(abs(sig.fdr$db.T1) >= dB),]
  print(paste("Control vs EoE.T1: ",nrow(sig.dB)))
}

hit_count(0.05, 0.01, dma.results)


###########################################
# DIFFERENTIAL METHYLATION ANALYSIS: DMRs #
###########################################

# Construct models & contrast matrices:

mod <- model.matrix(~0+dg.tpt+sex+age, data=meta)
corrfit <- duplicateCorrelation(mval, mod, block=meta$case.no)

contrastMatrix <- makeContrasts(
   T0 = dg.tptControl - dg.tptEoE.T0,
   T1 = dg.tptControl - dg.tptEoE.T1,
   eoe = dg.tptEoE.T0 - dg.tptEoE.T1,
   levels = mod)
                        
# Fit model & annotate with chromosome coordinates:
# NB, create EoE-T0 vs EoE-T1 model in order to double check that there are still no significant results for EoE-T0 vs EoE-T1 manually

T0.ann <- cpg.annotate("array", all.M, what="M", annotation=c(array="IlluminaHumanMethylationEPIC", annotation="ilm10b4.hg19"), analysis.type="differential", design=mod, contrasts=T, cont.matrix=contrastMatrix, fdr=0.01, coef="T0")
T1.ann <- cpg.annotate("array", all.M, what="M", annotation=c(array="IlluminaHumanMethylationEPIC", annotation="ilm10b4.hg19"), analysis.type="differential", design=mod, contrasts=T, cont.matrix=contrastMatrix, fdr=0.01, coef="T1")
eoe.ann <- cpg.annotate("array", all.M, what="M", annotation=c(array="IlluminaHumanMethylationEPIC", annotation="ilm10b4.hg19"), analysis.type="differential", design=mod, contrasts=T, cont.matrix=contrastMatrix, fdr=0.01, coef="eoe", block=meta$case.no, correlation=corrfit$consensus)

# Extract differentially methylated regions:
# NB, just for Controls vs EoE-T0 & vs EoE-T1 as know that there are no significant results for EoE-T0 vs EoE-T1:

T0.dmr <- dmrcate(T0.ann, lambda=1000, C=2)
T1.dmr <- dmrcate(T1.ann, lambda=1000, C=2)
#eoe.dmr <- dmrcate(eoe.ann, lambda=1000, C=2)

T0.dmr.results <- extractRanges(T0.dmr, genome="hg19")
T1.dmr.results <- extractRanges(T1.dmr, genome="hg19")
#eoe.dmr.results <- extractRanges(eoe.dmr, genome="hg19")
                        
# Save results:

readme.dma.3 <- data.frame(dataset=c("T0.ann","T1.ann","eoe.ann","T0.dmr","T1.dmr","T0.dmr.results","T1.dmr.results"),description=c("DMRcate ready annotated model for controls vs EoE.T0, excluding dups (n=19)","DMRcate ready annotated model for controls vs EoE.T1, excluding dups (n=18)","DMRcate ready annotated model for EoE.T0 vs EoE.T1, excluding dups (n=11)","DMRcate results for T0.ann","DMRcate results for T1.ann","hg19 annotated DMRcate results for T0.ann","hg19 annotated DMRcate results for T1.ann"))
save(T0.ann,T1.ann,eoe.ann,T0.dmr,T1.dmr,T0.dmr.results,T1.dmr.results,readme.dma.3,file="FINAL_ANALYSIS_RESULTS/DMA_RESULTS/EoE_DMRcate_DifferentialMethylationResults_Final.RData")

readme.dma.3

write.table(as.data.frame(T0.dmr.results), "FINAL_ANALYSIS_RESULTS/DMA_RESULTS/EoE_DMRcate_Control-vs-T0_DMRs_Final.txt", row.names=F, quote=F, sep="\t")
write.table(as.data.frame(T1.dmr.results), "FINAL_ANALYSIS_RESULTS/DMA_RESULTS/EoE_DMRcate_Control-vs-T1_DMRs_Final.txt", row.names=F, quote=F, sep="\t")

nrow(as.data.frame(T0.dmr.results))
nrow(as.data.frame(T1.dmr.results))


################
# SESSION INFO #
################

sessionInfo()
