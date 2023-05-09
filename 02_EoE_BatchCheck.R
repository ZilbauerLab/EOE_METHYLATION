# Check for significant batch effects
# 19Nov20 (fp215)
#
# Folders required:
#     EOE_METHYLATION/QC_PLOTS/04_BATCH_CORRECTION_CHECK
#
# Files required:
#     EOE_METHYLATION/OUTPUT_DATASETS/EoE_NormalisedFiltered_GRSet.RData (containing normalised, QCed dataset produced by 01_EoE_RGsetNormalisationQC.v2.6.R)


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
library(ggpubr)
library(RColorBrewer)
library(rafalib)
library(reshape)
library(scales)
library(sva)
library(BEclear)
library(data.table)

load("OUTPUT_DATASETS/EoE_EPICNormalisedData_GRSet.RData")
setwd("QC_PLOTS/04_BATCH_CORRECTION_CHECK")

readme.epic.grset

EPIC.GRset


######################################
# CREATE BECLEAR READY DATA AND INFO #
######################################

# Format key:

key <- data.frame(sample_id=samples$array.id,batch_id=samples$sentrix)
key2 <- data.table(key)

# Format dataset:

data <- data.table(feature=rownames(EPIC.beta), beta)
data <- melt(data=data, id.vars="feature", variable.name="sample", value.name="beta.value")
setkey(data, "feature", "sample")

# Calculate FDR p-values & median difference values for every probe in every batch:

pvals <- calcPvalues(data, key2, adjusted=T, method="fdr")
mdifs <- calcMedians(data, key2)

# Summarise the results:

summary <- calcSummary(median=mdifs, pvalues=pvals)


################################
# CALCULATE BATCH EFFECT SCORE #
################################

# No need for batch effect correction: BEscore < 0.02
# Defininitely needing batch correction: BEscore > 0.1
# Grey area (batch effect, but not severe): BEscore 0.02 - 0.1

score <- calcScore(beta, key2, summary, dir=getwd())

score

# Plot results by sample:

pdf("EoE_NormalisedFiltered_BEclearResults.pdf")
makeBoxplot(as.data.frame(beta), key, score, bySamples=T, col="standard", main="EoE Batch Effect Check", xlab="Batch", ylab="Beta value", scoreCol=T)
dev.off()


##########################
# SVA BATCH EFFECT CHECK #
##########################

# Check to see if top 2 calculated surrogate variables correlate with batch using sva:

mod <- model.matrix(~diagnosis,data=samples)
mod0 <- model.matrix(~1,data=samples)
svacheck <- sva(beta,mod,mod0,n.sv=2)

# Check to see if either surrogate variable correlates with batch (sentrix.id):

summary(lm(svacheck$sv ~ samples$sentrix))

# Plot results by batch:

samples$SV1 <- svacheck$sv[,1]
samples$SV2 <- svacheck$sv[,2]

require(gridExtra)

pdf("EoE_NormalisedFiltered_SVACheckResults.pdf")
p1 <- ggboxplot(samples, x="sentrix", y="SV1", width=0.3, color="sentrix", add="jitter", legend="none", title="Correlation between batch & SV1", ylab="SV1", xlab="Batch")
p2 <- ggboxplot(samples, x="sentrix", y="SV2", width=0.3, color="sentrix", add="jitter", legend="none", title="Correlation between batch & SV2", ylab="SV2", xlab="Batch")
grid.arrange(p1,p2,ncol=2)
dev.off()


################
# SESSION INFO #
################

sessionInfo()
