#################################################### Script for Creating STRUCTURE Plots  ########################################################

#using pophelper library as described in Francis 2016 - Molecular Ecology Resources
#manual found at: royfrancis.github.io/pophelper/#1_introduction
#STRUCTURE plots with all loci, loci only in HWE, outlier loci excluded, and only outlier loci

#################################################################################################################################################

######## Set-up ########

remove(list = ls())

#load libraries
library(devtools)
library(tidyverse)
library(gridExtra)
library(gtable)
#devtools::install_github('royfrancis/pophelper') #you will need to install this! (uncomment this line of code)
library(pophelper) #can't open until the previous line of code is run (one time only)

#read in data
SNPs_sfiles <- list.files(path = "STRUCTURE_Output/StructureResults/", full.names = TRUE)
SNPs_slist <- readQ(files = SNPs_sfiles, filetype = "structure") 

#create group labels 
grplab <- c(rep("Japan", 8), rep("Indonesia", 7), rep("Philippines", 10))#create group labels
meta.data <- data.frame(loc = grplab)
meta.data$loc <- as.character(meta.data$loc)

################################################################################################################################################

######## STRUCTURE Set-up ########

#summarize results
SNPs_table <- tabulateQ(qlist = SNPs_slist) #pulls parameters (including mean ln likelihood - mvli) from SNPs runs and puts into dataframe
SNPs_table_sum <- summariseQ(SNPs_table) #summarizes SNPs_table by K

#evanno method
SNPs_em <- evannoMethodStructure(data = SNPs_table_sum) #evanno method to pick best value of K
SNPs_em_plot <- evannoMethodStructure(data = SNPs_table_sum, exportplot = TRUE)

#clumpp
clumppExport(qlist = SNPs_slist, useexe = TRUE) #run clumpp to order clusters properly

######## Create STRUCTURE Plot ########

#read in CLUMPP data
SNPs_aligned_K2 <- readQ("Data/STRUCTURE_Output/CLUMPP/pop_K2-combined-aligned.txt")
SNPs_aligned_K3 <- readQ("Data/STRUCTURE_Output/CLUMPP/pop_K3-combined-aligned.txt")
SNPs_aligned_K4 <- readQ("Data/STRUCTURE_Output/CLUMPP/pop_K4-combined-aligned.txt")
SNPs_aligned_K5 <- readQ("Data/STRUCTURE_Output/CLUMPP/pop_K5-combined-aligned.txt")

#create plots
SNPs_K2 <- plotQ(SNPs_aligned_K2[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, 
                      clustercol = c("#2121D9", "#FF9329"),
                      showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
                      showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                      grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                      showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "all SNPS")

#repeat the SNPS_K2 line of code to create plots for K 3-5