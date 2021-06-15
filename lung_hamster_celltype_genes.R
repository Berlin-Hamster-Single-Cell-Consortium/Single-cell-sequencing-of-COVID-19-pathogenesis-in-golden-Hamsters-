library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(Seurat)
library(HDF5Array)
library(rhdf5)
library(hdf5r)
library(DelayedArray)
library(DelayedMatrixStats)
require(gridExtra)
library(forcats)
#see http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)
source("./summarySE.R")
library(lme4)

hamster.integrated <- readRDS("./ma_int.rds")




####################################
#Barplots and boxplots for genes in celltypes/virus levels etc. (Fig. 4 and 5)
#Barplots are with standard deviation and individual points
#Boxplots are for depicting expression values in those cells with at least one count for the gene
#Script also writes out individual values and p-values as text files


#Make barplots for positive for gene for several celltypes and all days with mean and sd and individual points
#This is for Fig. 4d and S5d
the_celltypes <- c("AT2", "MonocyticMacrophages")
#For the colors, use the original color as second color and color0 as the first
the_colors <- c(scales::seq_gradient_pal("#fdd8c6", "#f97e44", "Lab")(seq(0,1,length.out=5)), scales::seq_gradient_pal("#e9bdce", "#b7245c", "Lab")(seq(0,1,length.out=5)))

for (gene in c("Cxcl10", "Mx2", "Ccl2")) {
  df1 <- cbind.data.frame(FetchData(hamster.integrated, gene), hamster.integrated@meta.data$seurat_clusters, hamster.integrated@meta.data$celltype, str_replace_all(gene, "-", "_"), hamster.integrated@meta.data$SCoV2_load, hamster.integrated@meta.data$infection, hamster.integrated@meta.data$hamster)
  colnames(df1) <- c("expression", "cluster", "celltype", "gene", "load", "day", "hamster")
  df1 <- subset(df1, celltype %in% the_celltypes)
  
  a = df1 %>% group_by(day, celltype, hamster) %>% tally(name="tot") 
  b = df1 %>% group_by(day, celltype, hamster) %>% filter(expression>0) %>% tally(name="pos") 
  c =
    left_join(a , b , by = c('day', 'celltype', 'hamster')) %>% 
    replace(., is.na(.), 0) %>%
    mutate(fraction = pos / tot) %>%
    mutate(celltype = forcats::fct_relevel(celltype, the_celltypes))
  tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("day", "celltype"))
  write.table(c, paste("days", gene, "tsv", sep="."), sep="\t", quote=FALSE, row.names = FALSE)

  
  p1 <- ggplot(tgc, aes(x=celltype, y=fraction, fill=paste(celltype, day, sep="_")))+
    geom_bar(position=position_dodge(.8), stat="identity", width = 0.7)+
    geom_point(data=c, aes(x=celltype, y=fraction, fill=paste(celltype, day, sep="_")), position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size=0.2)+
    geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.8), size=0.25)+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    ggtitle(paste("fraction cells positive for", gene, ", mean and sd of three animals", sep=" "))+
    theme(legend.position = "none")+
    scale_fill_manual(values=the_colors)
  
  p2 <- ggplot()+
    geom_boxplot(data=subset(df1, expression>0), aes(x=celltype, y=expression, fill=paste(celltype, day, sep="_")), outlier.size=0.2, position=position_dodge(.8), width = 0.7, lwd=0.25, fatten=4)+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    ggtitle(paste("expression for", gene, "in cells positive for it",  sep=" "))+
    theme(legend.position = "none")+
    ylim(0, NA)+
    scale_fill_manual(values=the_colors)
  pdf(file=paste(gene, "percentpositive", "boxplotgt0", paste(the_celltypes, collapse="."), "pdf", sep="."), width=4, height=8, useDingbats=FALSE)
  grid.arrange(p1, p2, ncol=1)
  dev.off()
}

#Make barplots for positive for gene for several celltypes and all days with mean and sd, split in with/without virus
#Fig 5hi
the_celltypes <- c("AT2", "MonocyticMacrophages")
#For the colors, use the original color as second color and color0 as the first
the_colors <- c(scales::seq_gradient_pal("#fdd8c6", "#f97e44", "Lab")(seq(0,1,length.out=5)), scales::seq_gradient_pal("#e9bdce", "#b7245c", "Lab")(seq(0,1,length.out=5)))
for (gene in c("Map1lc3b", "Sqstm1")) {
  df1 <- cbind.data.frame(FetchData(hamster.integrated, gene), hamster.integrated@meta.data$seurat_clusters, hamster.integrated@meta.data$celltype, str_replace_all(gene, "-", "_"), hamster.integrated@meta.data$SCoV2_load, hamster.integrated@meta.data$infection, hamster.integrated@meta.data$hamster)
  colnames(df1) <- c("expression", "cluster", "celltype", "gene", "load", "day", "hamster")
  df1 <- subset(df1, celltype %in% the_celltypes)
  df1$virus <- ifelse(df1$load>0, "withvirus", "novirus")
  
  a = df1 %>% group_by(day, celltype, hamster, virus) %>% tally(name="tot") 
  b = df1 %>% group_by(day, celltype, hamster, virus) %>% filter(expression>0) %>% tally(name="pos") 
  c =
    left_join(a , b , by = c('day', 'celltype', 'hamster', 'virus')) %>% 
    replace(., is.na(.), 0) %>%
    mutate(fraction = pos / tot) %>%
    mutate(celltype = forcats::fct_relevel(celltype, the_celltypes))
  tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("day", "celltype", "virus"))
  
  expr <- list()
  #p-value calculation
  for (ct in the_celltypes) {
    for (d in c("d2", "d3", "d5")) {
      df2 <- subset(df1, day == d & celltype == ct)
      df2$expr <- ifelse(df2$expression>0, TRUE, FALSE)
      df2$celltype <- NULL
      df2$load <- NULL
      df2$expression <- NULL
      df2$cluster <- NULL
      df2$gene <- NULL
      df2$day <- NULL
      res <- try(p <- summary(glmer(df2$expr ~ df2$virus + (1 | df2$hamster), family=binomial))$coefficients[2,4])
      if(inherits(res, "try-error")) {
        p <- "NA"
      }
      expr[[paste(ct, d, sep="_")]] <- p
    }
  }
  pvalues <- as.data.frame(do.call(rbind, expr))
  write.table(pvalues, paste("pvalues_plusminusvirus_barplots", gene, "txt", sep="."), quote=FALSE, sep="\t", col.names=FALSE)

  expr <- list()
  #p-value calculation
  for (ct in the_celltypes) {
    for (d in c("d2", "d3", "d5")) {
      df2 <- subset(df1, day == d & celltype == ct)
       res <- try(p <- wilcox.test(x=subset(df2, virus == "novirus")$expression, y=subset(df2, virus == "withvirus")$expression)$p.value)
      if(inherits(res, "try-error")) {
        p <- "NA"
      }
      expr[[paste(ct, d, sep="_")]] <- p
    }
  }
  pvalues <- as.data.frame(do.call(rbind, expr))
  write.table(pvalues, paste("pvalues_plusminusvirus_boxplots", gene, "txt", sep="."), quote=FALSE, sep="\t", col.names=FALSE)
  write.table(c, paste("plusminusvirus", gene, "tsv", sep="."), sep="\t", quote=FALSE, row.names = FALSE)

  p1 <- ggplot(tgc, aes(x=celltype, y=fraction, fill=paste(celltype, day, sep="_"), color=virus))+
    geom_bar(position=position_dodge(.8), stat="identity", width = 0.7, size=0.25)+
    geom_point(data=c, aes(x=celltype, y=fraction, fill=paste(celltype, day, sep="_"), color=virus), position=position_jitterdodge(jitter.width = 0.04, dodge.width = 0.8), size=0.2)+
    geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.8), size=0.25)+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    ggtitle(paste("fraction positive for", gene, "with/without virus mean and sd of three animals", sep=" "))+
    scale_fill_manual(values=the_colors)+
    theme(legend.position = "none")+
    scale_color_manual(values=c("gray40", "gray20"))
  
  p2 <- ggplot()+
    geom_boxplot(data=subset(df1, expression>0), aes(x=factor(celltype, the_celltypes), y=expression, fill=paste(celltype, day, sep="_"), color=virus), outlier.size=0.2, position=position_dodge(.8), width = 0.7, lwd=0.25, fatten=4)+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    ggtitle(paste("expression for", gene, "in cells with/without positive for it",  sep=" "))+
    ylim(0, NA)+
    theme(legend.position = "none")+
    scale_fill_manual(values=the_colors)+
    scale_color_manual(values=c("gray40", "gray20"))
  pdf(file=paste(gene, "percentpositive_nowithvir", "boxplotgt0", paste(the_celltypes, collapse="."), "pdf", sep="."), width=4, height=6, useDingbats=FALSE)
  grid.arrange(p2, p1, ncol=1)
  dev.off()
}

#Fig 5j
#With error bars in the bottom and boxplot on top
quants <- c(0, 0.33, 0.66, 1)

#Colors: color4 is color from cell type assignement
#Go to https://www.color-hex.com/color/b7245c
#color3 is fourth color from left (first is original, which is color4), color0 is eigth
the_celltype <- "MonocyticMacrophages"
color0 = "#e9bdce"
color3 = "#cc658c"
color4 = "#b7245c"

the_celltype <- "AT2"
color0 = "#fdd8c6"
color3 = "#faa47c"
color4 = "#f97e44"

for (gene in c("Cxcl10")) {
  df1 <- cbind.data.frame(FetchData(hamster.integrated, gene), hamster.integrated@meta.data$seurat_clusters, hamster.integrated@meta.data$celltype, str_replace_all(gene, "-", "_"), hamster.integrated@meta.data$SCoV2_load, hamster.integrated@meta.data$infection, hamster.integrated@meta.data$hamster)
  colnames(df1) <- c("expression", "cluster", "celltype", "gene", "load", "day", "hamster")
  df1 <- subset(df1, celltype == the_celltype)
  df_d0 <- subset(df1, (day == "d0"))
  df_d0$level  <- 0
  expr_sub <- list()
  for (d in c("d2", "d3", "d5", "e14")) {
    df_temp <- subset(df1, (day == d))
    for (h in as.vector(unique(df_temp$hamster))) {
      df_temp1 <- subset(df_temp, (hamster == h & load > 0))
      df_temp1$level <- cut(df_temp1$load, breaks=quantile(df_temp1$load, probs=quants), label=FALSE, include.lowest = TRUE)
      df_temp1a <- subset(df_temp, (hamster == h & load == 0))
      df_temp1a$level <- 0
      expr_sub[[paste(d, h, sep="_")]] <- rbind.data.frame(df_temp1, df_temp1a)
    }
  }
  df_temp <- rbind.data.frame(df_d0, do.call(rbind.data.frame, expr_sub))
  
  a = df_temp %>% group_by(day, level, hamster) %>% tally(name="tot")
  b = df_temp %>% group_by(day, level, hamster) %>% filter(expression>0) %>% tally(name="pos")
  c = left_join(a , b , by = c('day', 'level', 'hamster')) %>% replace(., is.na(.), 0) %>% mutate(fraction = pos / tot)
  tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("day", "level"))
  
  expr <- list()
  #p-value calculation
  for (d in c("d2", "d3", "d5")) {
    for (l in unique(df_temp$level)) {
      for (ll in unique(df_temp$level)) {
        df2 <- subset(df_temp, day == d & (level == l | level == ll))
        df2$expr <- ifelse(df2$expression>0, TRUE, FALSE)
        df2$celltype <- NULL
        df2$load <- NULL
        df2$expression <- NULL
        df2$cluster <- NULL
        df2$gene <- NULL
        df2$day <- NULL
        res <- try(p <- summary(glmer(df2$expr ~ df2$level + (1 | df2$hamster), family=binomial))$coefficients[2,4])
        if(inherits(res, "try-error")) {
          p <- "NA"
        }
        expr[[paste(d, l, ll, sep="_")]] <- p
      }
    }
  }
  pvalues <- as.data.frame(do.call(rbind, expr))
  write.table(pvalues, paste("pvalues_levels_barplots", gene, "txt", sep="."), quote=FALSE, sep="\t", col.names=FALSE)
  write.table(c, paste("viruslevels", gene, "tsv", sep="."), sep="\t", quote=FALSE, row.names = FALSE)
 
  expr <- list()
  #p-value calculation
  for (d in c("d2", "d3", "d5")) {
    for (l in unique(df_temp$level)) {
      for (ll in unique(df_temp$level)) {
        df2 <- subset(df_temp, day == d & (level == l | level == ll))
        res <- try(p <- wilcox.test(x=subset(df2, level == l)$expression, y=subset(df2, level == ll)$expression)$p.value)
        if(inherits(res, "try-error")) {
          p <- "NA"
        }
        expr[[paste(d, l, ll, sep="_")]] <- p
      }
    }
  }
  pvalues <- as.data.frame(do.call(rbind, expr))
  write.table(pvalues, paste("pvalues_levels_boxplots", gene, "txt", sep="."), quote=FALSE, sep="\t", col.names=FALSE)
  
  
   p1 <- ggplot(tgc, aes(x=day, y=fraction, fill=as.factor(level)))+
    geom_bar(position=position_dodge(.8), stat="identity", width = 0.7, size=0.25)+
    geom_point(data=c, aes(x=day, y=fraction, fill=as.factor(level)), position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size=0.2)+
    geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.8), size=0.25)+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    theme(legend.position = "none")+
    ggtitle(paste("fraction positive", gene, "in virus levels mean and sd of three animals", sep=" "))+
    scale_fill_manual(values=c(color0, scales::seq_gradient_pal(color3, color4, "Lab")(seq(0,1,length.out=length(quants)-1))))

   df_d0 <- subset(df1, (day == "d0"))
   df_d0$level  <- 0
   expr_sub <- list()
   for (d in c("d2", "d3", "d5", "e14")) {
     df_temp <- subset(df1, (day == d))
     df_temp1 <- subset(df_temp, load > 0)
     df_temp1$level <- cut(df_temp1$load, breaks=quantile(df_temp1$load, probs=quants), label=FALSE, include.lowest = TRUE)
     df_temp1a <- subset(df_temp, load == 0)
     df_temp1a$level <- 0
     expr_sub[[d]] <- rbind.data.frame(df_temp1, df_temp1a)
   }
   df_temp <- rbind.data.frame(df_d0, do.call(rbind.data.frame, expr_sub))
   p2 <- ggplot()+
     geom_boxplot(data=subset(df_temp, expression>0), aes(x=day, y=expression, fill=as.factor(level)), outlier.size=0.2, position=position_dodge(.8), width = 0.7, lwd=0.25, fatten=4)+
     theme_bw()+
     theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
     theme(legend.position = "none")+
     ggtitle(paste("expression for", gene, "in", the_celltype, "cells in virus levels positive for it", sep=" "))+
     ylim(0, NA)+
     scale_fill_manual(values=c(color0, scales::seq_gradient_pal(color3, color4, "Lab")(seq(0,1,length.out=length(quants)-1))))
   pdf(file=paste(gene, "percentpositive", "boxplotgt0", the_celltype, "pdf", sep="."), width=6, height=8, useDingbats=FALSE)
   grid.arrange(p1, p2, ncol=1)
   dev.off()
}
