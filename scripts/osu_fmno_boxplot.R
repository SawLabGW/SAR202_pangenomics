#/usr/bin/env Rscript
## Script used to generate Figure 2B
## example
# cd ~/cgrb_workspace/SAR202/refined_final_set/stats
# Rscript ~/workspace/repositories/osu_work/scripts/osu_fmno_boxplot.R -c enzyme_abundances.txt -p info -o f1.pdf

## Load libraries
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("yarrr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggthemes"))

setwd("~/cgrb_workspace/SAR202/refined_final_set/stats")


#tib <- read_tsv("enzyme_abundances_original_bin_names.txt")
tib <- read_tsv("enzyme_abundances_with_AslA.roman.txt")
df <- as.data.frame(tib)

## remove "Group" to make names shorter
#df[,2] <- as.character(gsub("Group", "", as.matrix(df[,2])))

#pdf(file = args$output,   # The directory you want to save the file in
#    width = 11, # The width of the plot in inches
#    height = 8)

## set variables for title and palette
#t <- args$title
#p <- args$palette

# pirateplot(formula = abundance ~ group + enzyme,
#            data = df,
#            main = "test",
#            pal = p,
#            theme = 2,
#            avg.line.col = "black", # avg line col
#            avg.line.o = 1, # Average line
#            bar.f.col = gray(.8), # bar filling color
#            bean.b.o = .2,
#            bean.f.o = .6, # Bean fill
#            inf.f.o = .7, # Inference fill
#            inf.b.o = .8, # Inference border
#            #inf.disp = "bean",
#            #inf.f.col = "white", # Inf fill col
#            inf.b.col = "black", # Inf border col
#            jitter.val = 0.06,
#            point.o = .3, # Points
#            point.bg = "black",
#            point.col = "black",
#            point.cex = .8,
#            point.pch = 21,
#            quant = c(.1, .9), # 10th and 90th quantiles
#            quant.col = "black",
#            sortx = "alphabetical",
#            #decreasing = TRUE,
#            ylab = "% of total genes",
#            xlab = "SAR202 subgroup")

p <- ggplot(data = df, aes(x=group, y=abundance))
p <- p + geom_boxplot(aes(fill=enzyme),
                      outlier.size=0.5,
                      outlier.colour="#6A5ACD")

p <- p + scale_fill_manual(values=c("#EEB4B4", "#D1EEEE", "#CDC673", "#EED2EE", "#7CCD7C"),
                           labels=c("arylsulfatases", "dehydrogenases","enolases","FMNOs",
                                    "RHDs"))
p <- p + geom_jitter(width = 0.25, alpha=0.5, size=0.5, colour="#7F7F7F")
p <- p + facet_wrap( ~ enzyme, scales="free", ncol=3)
p <- p + xlab("SAR202 subgroups") + ylab("% abundance of total genes") #+ ggtitle("SAR202 enzyme abundances")
p <- p + guides(fill=guide_legend(title="Enzymes"))

#p + theme_light()
#p + theme_tufte() 
#p + theme_hc()
p + theme_minimal() + 
  theme( axis.text = element_text( size = 12 ),
         axis.text.x = element_text( size = 12 ),
         axis.title = element_text( size = 12, face = "bold" ),
         legend.position=c(0.85,0.2), legend.text = element_text(size = 12),
         strip.text = element_text(size = 12))
