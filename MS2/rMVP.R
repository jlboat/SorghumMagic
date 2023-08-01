library(rMVP)
setwd("/zfs/tillers/nkumar2/GWAS_Yeild_Traits")
df <- read.csv("phenotypes.csv", header=T)
write.table(df, file=" phenotypes.txt", quote=F, row.names = F)
MVP.Data(fileVCF="MAGIC.quality.renamed.vcf",
         filePhe="phenotypes.csv",
         sep.phe = ",",
         fileKin=F,
         filePC=F,
         out="mvp.vcf"
)
genotyp_dat <- attach.big.matrix("mvp.vcf.geno.desc")
phenotypic_dat <- read.table("mvp.vcf.phe",head=TRUE)
map_info <- read.table("mvp.vcf.geno.map" , head = TRUE)
GWAS_imMVP <- MVP(
  phe = phenotypic_dat,
  geno=genotyp_dat,
  map=map_info,
  method=c("MLM", "FarmCPU"),
  nPC.MLM=3,
  threshold=0.05,
)

