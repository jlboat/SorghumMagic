source("http://zzlab.net/GAPIT/gapit_functions.txt")
library("scatterplot3d")
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(EMMREML)
library(compiler)

setwd <- ("~/Desktop/GWAS1")
myY <- read.csv("~/Desktop/GWAS1/phenotypes.csv", header = TRUE)
myG <-read.table(file="~/Desktop/GWAS1/MAGIC_DArT_v3.filtered.imputed.biallelic.hmp.txt", header = F, sep="\t", comment.char = "%")

#first column is individual ID, the 'i' column is the phenotype of interest
for (i in seq(2,5,1)){
    myGAPIT_MLM <- GAPIT(
                         Y=myY[,c(1,i)],
                         G=myG,
                         PCA.total=3,
                         model=c("MLM"),
                         ncpus=2)
    pheno <- colnames(myY)[i]
    file.rename("GAPIT.Filter_GWAS1_results.txt", 
                paste0("GAPIT.", 
                       pheno,
                       ".Filter_GWAS1_results.txt"))
}












