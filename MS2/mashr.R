rm(list=ls())
#install.packages("devtools")
#devtools::install_github("stephenslab/mashr@v0.2-11")
library(foreach)
library(iterators)
library(parallel)
library(ashr)
library(mashr)
library(mixsqp)
#library(corrplot)
#library(ggplot2)
library(data.table)
library(doParallel)
registerDoParallel(cores = 16)
library(RhpcBLASctl)
blas_set_num_threads(16)
mc.core = 16

setwd("/zfs/tillers/nkumar2/rMVP-GWAS_Yield_Traits")

# Merge specific columns into new files with first column as standard or by the first column
# Here input is GWAS results from rMVP. Results from other softwares can be used, make sure to correctly assign the column numbers of respective beta/effect estimates and SE
# May change this line to run FarmCPU.csv, etc.
ls <- list.files(pattern = ".MLM.csv")
x <- fread(paste(ls[1], sep = "")) 
# x[,"SNP"] <- paste(x$chr, x$ps, sep="_")
datBeta <- x[,c(1, 6)]
datSE <- x[,c(1, 7)]
colnames(datBeta) <- c("SNP", ls[1])
colnames(datSE) <- c("SNP", ls[1])

datPtopMLM <- data.frame()
for (i in 2:length(ls)){          #i in 1:176 or i in 1:length(ls) # make sure there are no other csv file in the folder apart from GWAS results that you want to include in this analysis
  print(i)
  x <- fread(paste(ls[i], sep = ""))
  x <- na.omit(x)
  # x[,"SNP"] <- paste(x$chr, x$ps, sep="_")
  datBeta <- merge(datBeta, x[,c(1, 6)], by="SNP")
  colnames(datBeta)[i+1] <- ls[i]
  datSE <- merge(datSE, x[,c(1, 7)], by = "SNP")
  colnames(datSE)[i+1] <- ls[i]
  pmlm <- x[order(x[,1])[1], c(1,2,3,4,5,8)] # chr ps allele1 allele0 12?
  pmlm$trait <- ls[i]
  colnames(pmlm)[6] <- "p"
  datPtopMLM <- rbind(datPtopMLM, pmlm)
  print(paste(ls[i], "Done!"))
}
saveRDS(datPtopMLM, file="all_top_MLM.rds")
 
datBeta <- na.omit(datBeta)
datSE <- na.omit(datSE)
# 
# if you want to save your datBeta and datSE files so that you dont have to run above command again
saveRDS(datBeta, "datBeta.rds", compress = T)
saveRDS(datSE, "datSE.rds", compress = T)
###### load files #
# use below command and start from this step after loading packages if you already have datBeta and datSE files

#datBeta <- readRDS("datBeta.rds")
#datSE <- readRDS("datSE.rds")
# MashR does NOT accept dataframe so change dataframe to matrix format
#rownames(datBeta) <- datBeta$SNP
#rownames(datSE) <- datSE$SNP
row_beta <- datBeta$SNP
row_se <- datSE$SNP

datBeta <- datBeta[,-c(1)]
datSE <- datSE[,-c(1)]

#na <-  -c(which(is.na(datBeta)))
#datBeta <- datBeta[na,]
#datSE <- datSE[na,]

datBeta <- as.matrix(datBeta)
datSE <- as.matrix(datSE)

rownames(datBeta) <- row_beta
rownames(datSE) <- row_se
################# MASH #####################  
# identify a subset of strong tests
m.1by1 = mash_1by1(mash_set_data(datBeta,datSE))
strong.subset = get_significant_results(m.1by1, 0.1)

saveRDS(strong.subset, "strong.subset.rds", compress = T)
###### load files #
# strong.subset <- readRDS("strong.subset.rds")
# identify a random subset of 90000 tests
random.subset = sample(1:nrow(datBeta), 1000)
saveRDS(random.subset, "random.subset.rds", compress = T)
#random.subset <- readRDS("./scr/176random.subset.rds")
# Correlation structure
data.temp = mash_set_data(datBeta[random.subset,],datSE[random.subset,])
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)
data.random = mash_set_data(datBeta[random.subset,],datSE[random.subset,],V=Vhat)
data.strong = mash_set_data(datBeta[strong.subset,],datSE[strong.subset,], V=Vhat)
saveRDS(data.random, "data.random.rds", compress = T)
saveRDS(data.strong, "data.strong.rds", compress = T)
###### load files #
# data.random<- readRDS("data.random.rds")
# data.strong <- readRDS("data.strong.rds")

# Data driven covariances
U.pca = cov_pca(data.strong, 5)
U.ed = cov_ed(data.strong, U.pca)
# Fit mash model (estimate mixture proportions)
U.c = cov_canonical(data.random) # remove if do not want to add canonical and run faster
m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1) # remove if do not want to add canonical and run faster

# Fit mash model (estimate mixture proportions)
#m = mash(data.random, Ulist = c(U.ed), outputlevel = 1)

saveRDS(m, "m.rds", compress = T)

# Compute posterior summaries
m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)

# Sharing
print(get_pairwise_sharing(m2))
#print(get_pairwise_sharing(m.c, factor=0))
saveRDS(m2, "m2.rds", compress = T)
sessionInfo()

