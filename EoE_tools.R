# PLOTTING COLOUR SCHEMES, THEMES & FUNCTIONS
# heat.scree function adapted from Rachel Edger's Heat_scree_scree_plot_generic.R (https://github.com/redgar598, 01Apr19)
# 23Nov20 (fp215)

library(ggplot2)
library(RColorBrewer)
library(scales)
library(gridExtra)
library(reshape)
library(nlme)


# COLOUR SCHEMES

myColours_diagnosis_eoe <- c("grey75","darkred")
colour_possibilities_diagnosis_eoe <- c("Control","EoE")
names(myColours_diagnosis_eoe) <- colour_possibilities_diagnosis_eoe
fillscale_diagnosis_eoe <- scale_fill_manual(name="Diagnosis", values = myColours_diagnosis_eoe, drop = T)
colscale_diagnosis_eoe <- scale_colour_manual(name="Diagnosis", values = myColours_diagnosis_eoe, drop = T)

myColours_timept_eoe <- c("black","dodgerblue","dodgerblue4")
colour_possibilities_timept_eoe <- c("T0","T1","T2")
names(myColours_timept_eoe) <- colour_possibilities_timept_eoe
fillscale_time_eoe <- scale_fill_manual(name="Time Point", values=myColours_timept_eoe, drop=T)
colscale_time_eoe <- scale_colour_manual(name="Time Point", values=myColours_timept_eoe, drop=T)

myColours_eoe <- c("grey75","darkblue","darkorchid","darkred")
colour_possibilities_eoe <- c("Control","EoE.T0","EoE.T1","EoE.T2")
names(myColours_eoe) <- colour_possibilities_eoe
fillscale_eoe <- scale_fill_manual(name="DiagnosisTimepoint",values=myColours_eoe,drop=T)
colscale_eoe <- scale_colour_manual(name="DiagnosisTimepoint",values=myColours_eoe,drop=T)

myColours_heatscree <- c("#084594","#4292c6","#9ecae1","#deebf7")
colour_possibilities_heatscree <- c("<=0.001","<=0.01","<=0.05",">0.05")
names(myColours_heatscree) <- colour_possibilities_heatscree
fillscale_heatscree <- scale_fill_manual(name="P value", values = myColours_heatscree, drop=T)
colscale_heatscree <- scale_colour_manual(name="P value", values = myColours_heatscree, drop=T)


# THEMES

myTheme <- theme(axis.text=element_text(size=12), axis.title=element_text(size=14), strip.text.x = element_text(size=12), legend.text=element_text(size=12), legend.title=element_text(size=14))


# HEAT SCREE PLOT FUNCTION

heat.scree <- function(Loadings, Importance, right_marg, left_marg){

    if(missing(right_marg)){
        right_marg=2.25}

    if(missing(left_marg)){
        left_marg=1}

    pca.df <- data.frame(variance=Importance, PC=seq(1:length(Importance)))
    
    scree <- ggplot(pca.df[which(pca.df$PC<=(PCs_to_view)),],aes(PC,variance))+
        geom_bar(stat="identity",color="black",fill="grey")+
        theme_bw()+
        theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title=element_text(size=15),plot.margin=unit(c(1.25,2,-0.2,3),"cm"))+
        ylab("Variance")+
        scale_x_continuous(breaks=seq(1,PCs_to_view,1))+
        xlab("")

# Correlate meta data with PCs:

    aov.PC.meta <- lapply(1:ncol(meta_categorical),function(covar) sapply(1:ncol(Loadings), function(PC) summary(aov(Loadings[,PC] ~ meta_categorical[,covar]))[[1]]$"Pr(>F)"[1]))
    cor.PC.meta <- lapply(1:ncol(meta_continuous), function(covar) sapply(1:ncol(Loadings), function(PC) (cor.test(Loadings[,PC], as.numeric(meta_continuous[,covar]), alternative="two.sided", method="kendall", na.action=na.omit)$p.value)))

    names(aov.PC.meta) <- colnames(meta_categorical)
    names(cor.PC.meta) <- colnames(meta_continuous)
    aov.PC.meta <- do.call(rbind, aov.PC.meta)
    cor.PC.meta <- do.call(rbind, cor.PC.meta)
    aov.PC.meta <- rbind(aov.PC.meta, cor.PC.meta)
    aov.PC.meta <- as.data.frame(aov.PC.meta)

# Adjust and reshape:

    aov.PC.meta.adjust <- aov.PC.meta[,1:ncol(aov.PC.meta)]
    avo <- aov.PC.meta.adjust[,1:PCs_to_view]
    avo.heat.no <- apply(avo,2,as.numeric)
    avo.heat <- as.data.frame(avo.heat.no)
    avo.heat$meta <- rownames(avo)
    avo.heat.melt <- melt(avo.heat, id=c("meta"))
    meta.var.order <- unique(avo.heat.melt$meta)[rev(ord)]
    avo.heat.melt$meta <- factor(avo.heat.melt$meta, levels=meta.var.order)

# Highlight if significant:

    avo.heat.melt$Pvalue <- sapply(1:nrow(avo.heat.melt), function(x) if(avo.heat.melt$value[x]<=0.001){"<=0.001"}else{
        if(avo.heat.melt$value[x]<=0.01){"<=0.01"}else{
            if(avo.heat.melt$value[x]<=0.05){"<=0.05"}else{">0.05"}}})

    levels(avo.heat.melt$variable) <- sapply(1:PCs_to_view, function(x) paste("PC",x,sep=""))

# Create plot:

    heat <- ggplot(avo.heat.melt, aes(variable, meta, fill=Pvalue))+
        geom_tile(color="black", size=0.5)+
        theme_gray(8)+fillscale_heatscree+
        theme(axis.text=element_text(size=10,color="black"),axis.text.x=element_text(),axis.title=element_text(size=15),legend.text=element_text(size=10),legend.title=element_text(size=10),legend.position=c(1.15,0.75),legend.justification=c(1,1),plot.margin=unit(c(-0.3,right_marg,1,left_marg),"cm"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
        xlab("Principal Component")+
        ylab(NULL)

    grid.arrange(scree, heat, ncol=1,heights = c(3, 4))
}

