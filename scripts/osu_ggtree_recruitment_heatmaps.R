#/usr/bin/env Rscript

## Load libraries
library("tidyverse")
library("ggplot2")
library("ggtree")

setwd("~/cgrb_workspace/SAR202/refined_final_set/recruitment/")

## load tree
tree <- read.newick("chains14_b300_10_2600.con.pruned.tre")

## load group info
info <- read.csv("new_groups_110219.csv")

## set colors for groups
cols <- c(Ia='#EED2EE', Ib='#CDB6CD', Ic='#CD96CD', II='#607B8B',
          IIIa='#EECFA1', IIIb='#CDB38B', IIIc='#8B795E', IV='#8FBC8F', 
          V='#FA8072', VI='#008080', VII='#C70039')

## plot recruitment data for trenches
trenches <- read_tsv("mod_rec_trenches.txt")
df <- as.data.frame(trenches)
rownames(df) <- df$bins
df$bins <- NULL

t1 <- ggtree(tree, color="black", branch.length="none") %<+% info +
  theme_tree() +
  geom_tippoint(aes(color=group, shape=15)) + 
  scale_shape_identity() +
  scale_color_manual(values=cols) + coord_flip() + scale_x_reverse()

t2 <- gheatmap(t1, 
               df, 
               width=10, 
               low="white", 
               high="black",
               #color="#EBEBEB",
               colnames_position="top",
               font.size=3,
               offset=0,
               colnames_angle=0,
               colnames_offset_y=1,
               hjust=0)

plot(t2)

## plot recruitment data for TARA
tara <- read_tsv("rec_tara.txt")
df <- as.data.frame(tara)
rownames(df) <- df$bins
df$bins <- NULL

tara1 <- ggtree(tree, color="black", branch.length="none") %<+% info +
  theme_tree() +
  geom_tippoint(aes(color=group, shape=15)) + 
  scale_shape_identity() +
  scale_color_manual(values=cols) + coord_flip() + scale_x_reverse()

tara2 <- gheatmap(tara1, 
               df, 
               width=10, 
               low="white", 
               high="black",
               #color="#EBEBEB",
               colnames_position="top",
               font.size=2.5,
               offset=0,
               colnames_angle=0,
               colnames_offset_y=1,
               hjust=0)

plot(tara2)

## plot recruitment data for BATS
bats <- read_tsv("rec_bats.txt")
df <- as.data.frame(bats)
rownames(df) <- df$bins
df$bins <- NULL

bats1 <- ggtree(tree, color="black", branch.length="none") %<+% info +
  theme_tree() +
  geom_tippoint(aes(color=group, shape=15)) + 
  scale_shape_identity() +
  scale_color_manual(values=cols) + coord_flip() + scale_x_reverse()

bats2 <- gheatmap(bats1, 
                  df, 
                  width=10, 
                  low="white", 
                  high="black",
                  #color="#EBEBEB",
                  colnames_position="top",
                  font.size=4,
                  offset=0,
                  colnames_angle=0,
                  colnames_offset_y=1,
                  hjust=0)

plot(bats2)
