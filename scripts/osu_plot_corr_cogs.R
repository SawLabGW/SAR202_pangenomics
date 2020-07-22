library(tidyverse)
library(corrplot)
library(RColorBrewer)

## Script used to make Figure 3

## make sure to mount folder from server before doing this
setwd("~/cgrb_workspace/SAR202/refined_final_set/rpsblast/results")

## load the data frame produced by osu_plot_cog_heatmap.py script
tib0 <- read_tsv("data_frame.tab")
tib1 <- read_tsv("group1.df")
tib2 <- read_tsv("group2.df")
tib3 <- read_tsv("group3.df")
tib4 <- read_tsv("group4.df")
tib5 <- read_tsv("group5.df")
#tib24 <- read_tsv("test.df")

df0 <- as.data.frame(tib0)
df1 <- as.data.frame(tib1)
df2 <- as.data.frame(tib2)
df3 <- as.data.frame(tib3)
df4 <- as.data.frame(tib4)
df5 <- as.data.frame(tib5)
#df24 <- as.data.frame(tib24)

rownames(df0) <- df0$bins
df0$bins <- NULL

rownames(df1) <- df1$bins
df1$bins <- NULL

rownames(df2) <- df2$bins
df2$bins <- NULL

rownames(df3) <- df3$bins
df3$bins <- NULL

rownames(df4) <- df4$bins
df4$bins <- NULL

rownames(df5) <- df5$bins
df5$bins <- NULL

#rownames(df24) <- df24$bins
#df24$bins <- NULL

#col <- colorRampPalette(c("#548B54", "#7CCD7C", "#FFFFFF", "#EE9988", "#BB4444"))
#col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#7CCD7C", "#548B54"))
col <- colorRampPalette(c("#8B814C", "#CDBE70", "#FFFFFF", "#5CACEE", "#36648B"))

M0 <- cor(df0)
pdf("R_corrplot_all.pdf")
#corrplot(M0, order="hclust", col=col(n=10), tl.col="#808080", tl.cex=.6)
#corrplot(M0, order="hclust", hclust.method="ward", col=col(n=10), tl.col="#808080", tl.cex=.6)
#corrplot(M0, order="hclust", hclust.method="ward.D", col=col(n=10), tl.col="#808080", tl.cex=.6)
corrplot(M0, 
         order="hclust", 
         hclust.method="complete", 
         col=col(n=10), 
         tl.col="#808080", 
         tl.cex=.6,
         insig="blank")
dev.off()

M1 <- cor(df1)
pdf("R_corrplot_group1.pdf")
#corrplot(M1, order="hclust", col=col(n=10), tl.col="#808080", tl.cex=.6)
corrplot(M1, order="hclust", hclust.method="complete", col=col(n=10), tl.col="#808080", tl.cex=.6)
dev.off()

M2 <- cor(df2)
pdf("R_corrplot_group2.pdf")
#corrplot(M2, order="hclust", col=col(n=10), tl.col="#808080", tl.cex=.6)
corrplot(M2, order="hclust", hclust.method="complete", col=col(n=10), tl.col="#808080", tl.cex=.6)
dev.off()

M3 <- cor(df3)
pdf("R_corrplot_group3.pdf")
#corrplot(M3, order="hclust", col=col(n=10), tl.col="#808080", tl.cex=.6)
corrplot(M3, order="hclust", hclust.method="complete", col=col(n=10), tl.col="#808080", tl.cex=.6)
dev.off()

M4 <- cor(df4)
pdf("R_corrplot_group4.pdf")
#corrplot(M4, order="hclust", col=col(n=10), tl.col="#808080", tl.cex=.6)
corrplot(M4, order="hclust", hclust.method="complete", col=col(n=10), tl.col="#808080", tl.cex=.6)
dev.off()

M5 <- cor(df5)
pdf("R_corrplot_group5.pdf")
#corrplot(M5, order="hclust", col=col(n=10), tl.col="#808080", tl.cex=.6)
corrplot(M5, order="hclust", hclust.method="complete", col=col(n=10), tl.col="#808080", tl.cex=.6)
dev.off()

#corrplot(M, order="hclust", col=brewer.pal(n=10, name="PuOr"), tl.col="#808080", tl.cex=.6)
#corrplot(M, order="hclust", type="upper", col=brewer.pal(n=10, name="PuOr"), tl.col="#808080", tl.cex=.6)


