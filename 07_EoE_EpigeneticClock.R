# Epigenetic clock analysis
# 25Nov20 (fp215)
#
# Folders required:
#     EOE_METHYLATION/EPIGENETIC_CLOCK
#
# Files required:
#     OUTPUT_DATASETS/EoE_RawData_RGSet.RData (raw data prior to normalisation and filtering)
#     EOE_METHYLATION/EPIGENETIC_CLOCK/Horvath_datMiniAnnotation.csv (Horvath CpG file: "datMiniAnnotation3.csv")
#
# DNA methylation age of human tissues and cell types, Horvath, Genome Biol 14(10):R115 (2013)
#     https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4015143/pdf/gb-2013-14-10-r115.pdf
# Systematic evaluation of DNA methylation age estimation with common preprocessing methods and the Infinium MethylationEPIC BeadChip array, McEwan et al, Clin Epigenetics 10(1):123 (2018)
#     https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-018-0556-2


#######################
# WORKING ENVIRONMENT #
#######################

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
library(lme4)

setwd("~/EOE_METHYLATION")
load("OUTPUT_DATASETS/EoE_RawData_RGSet.RData")

readme.rgset


##################################################
# QC, NORMALISE AND FILTER INCLUDING CHR X AND Y #
##################################################

pval <- detectionP(RGset.raw)
pval.mean <- colMeans(detectionP(RGset.raw))
GRset <- preprocessFunnorm(RGset, nPCs=2, bgCorr=T, dyeCorr=T)
beta <- getBeta(GRset)

na.count <- sapply(1:nrow(beta), function(y) length(which(is.na(beta[y,]))))
na.count <- which(na.count < (ncol(beta)*0.05))
beta <- beta[na.count,]

fails <- pval > 0.05
fail.p <- names(which(rowMeans(fails) > 0.01))
fail.id <- names(which(colMeans(fails) > 0.01))
beta <- beta[which(!(rownames(beta) %in% fail.p)),]
beta <- beta[,which(!(colnames(beta) %in% fail.id))]

# Check beta values and sample info datasets match up

all(colnames(beta)==rownames(samples$assay.id))


###############################################
# CREATE HORVATH EPIGENETIC CLOCK READY INPUT #
###############################################

# NB, removing outliers (SZ1.T0 and SZ4.T2)

info <- data.frame(Sample.ID=samples$sample.id,Age=round(samples$age,1),Tissue=rep("Head+Neck",31),Female=ifelse(samples$sex=="F",0),Diagnosis=samples$dg.tpt,eos=samples$eos,Case.ID=samples$case.no)
info <- subset(info, Sample.ID!="SZ1.T0" & Sample.ID!="SZ4.T2")
write.csv(info, file="EPIGENETIC_CLOCK/EoE_ClockInput_AllSampleInfo_n29.csv", row.names=F, quote=F)

data <- as.data.frame(beta)
colnames(data) <- samples$sample.id
data$ProbeID <- rownames(data)
data <- data[,c(32,1:31)]
data$SZ1.T0 <- NULL
data$SZ4.T2 <- NULL
write.csv(data, file="EoE_ClockInput_AllSampleBetaValues_n29.csv", row.names=F, quote=F)

# Subset to just the CpGs used by the Horvath clock

cpgs <- read.table("EPIGENETIC_CLOCK/Horvath_datMiniAnnotation.csv", header=T, sep=",", stringsAsFactors=F)
list <- cpgs$Name
list <- data.frame(ProbeID=list)
horvath <- join(list, data, type="left", match="all")
write.csv(horvath, file="EPIGENETIC_CLOCK/EoE_AllSampleHorvathBetaValues_n29.csv", row.names=F, quote=F)


############################
# HORVATH EPIGENETIC CLOCK #
############################

# Run from online tool: https://dnamage.genetics.ucla.edu/new
# Use "info" and "horvath" input created above, selecting the "Normalize Data" option, following guidelines


###################
# ANALYSE RESULTS #
###################

# Simplified table of results produced from online clock: EPIGENETIC_CLOCK/EoE_AllSamples_HorvathResults_n29.txt

results <- read.table("EPIGENETIC_CLOCK/EoE_AllSamples_HorvathResults_n29.txt", header=T, sep="\t", stringsAsFactors=F)
results <- results[c("Sample.ID","Case.ID","Age","DNAmAge","AgeAccelerationDiff","AgeAccelerationResidual","Female","Diagnosis","eos","meanMethBySample","minMethBySample","maxMethBySample","corSampleVSgoldstandard","meanAbsDifferenceSampleVSgoldstandard","meanXchromosome","ProbabilityFrom.Head.Neck","ProbabilityFrom.Buccal")]
results <- rename(results, c("AgeAccelerationDiff"="Diff","AgeAccelerationResidual"="Residual","meanMethBySample"="meanMeth","minMethBySample"="minMeth","maxMethBySample"="maxMeth","corSampleVSgoldstandard"="corSample","meanAbsDifferenceSampleVSgoldstandard"="meanAbsDiffSample","meanXchromosome"="meanX","ProbabilityFrom.Head.Neck"="PrHeadNeck","ProbabilityFrom.Buccal"="PrBuccal"))

# Calculate & plot overall correlation between known chronological age & epigenetic predicted age in all samples

age.cor <- cor.test(results$Age,results$DNAmAge)

age.cor

pdf("EPIGENETIC_CLOCK/EoE_HorvathClock_ChronologicalAge-v-EpigeneticAge_DiagnosisTimepoint.pdf")
ggplot(results, aes(Age, DNAmAge, fill=Diagnosis))+geom_point(size=2,color="black",shape=21)+xlim(0,30)+ylim(0,30)+theme_bw()+scale_fill_manual("Diagnosis",values=c("grey75","plum","purple"))+geom_abline(slope=1, intercept=0, color="grey")+annotate("text",x=25,y=5,label=paste0("r = ",round(age.cor$estimate,2)",\np = ",round(age.cor$p.value,3)),size=5)+xlab("Chronological Age")+ylab("DNAm Age")
dev.off()

pdf("EPIGENETIC_CLOCK/EoE_HorvathClock_ChronologicalAge-v-EpigeneticAge_Eosinophils-per-hpf.pdf")
ggplot(results, aes(Age, DNAmAge, fill=eos))+geom_point(size=2,color="black",shape=21)+xlim(0,30)+ylim(0,30)+theme_bw()+scale_fill_manual("Eosinophils",values=c("gold","darkorange","red","grey95"),labels=c("15-20 per hpf","20-50 per hpf","Over 50 per hpf","Under 15 per hpf"))+geom_abline(slope=1, intercept=0, color="grey")+annotate("text",x=25,y=5,label=paste0("r = ",round(age.cor$estimate,2),"\np = ",round(age.cor$p.value,3)),size=5)+xlab("Chronological Age")+ylab("DNAm Age")
dev.off()

pdf("EPIGENETIC_CLOCK/EoE_HorvathClock_Diagnosis-v-AgeAccelerationResidual.pdf")
ggplot(results, aes(Diagnosis, Residual, fill=Diagnosis))+geom_boxplot(color="black",outlier.shape=NA)+geom_point(shape=21,position=position_jitter(w=0.15))+theme_bw()+scale_fill_manual("Diagnosis",values=c("grey75","plum","purple"))+ylab("Age Acceleration Residual")
dev.off()

# Compare sample groups

T0.results <- subset(results, Diagnosis!="EoE.T1")
T1.results <- subset(results, Diagnosis!="EoE.T0")
eoe.results <- subset(results, Diagnosis!="Control")

summary(aov(results$Residual~results$Diagnosis))

summary(aov(T0.results$Residual~T0.results$Diagnosis))

summary(aov(T1.results$Residual~T1.results$Diagnosis))

eoe.results <- subset(results, Diagnosis!="Control")

mod1 <- lmer(T0.results$Residual~T0.results$Diagnosis+(1|T0.results$Case.ID))
mod2 <- lmer(T0.results$Residual~(1|T0.results$Case.ID))
T0.lrt <- anova(mod1,mod2)

T0.lrt
