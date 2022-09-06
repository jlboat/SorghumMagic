library(qtl2)
library(LinkageMapView)

magic <- ("~/Desktop/MAGIC1")
magic <- read_cross2("~/Desktop/MAGIC1/magic.yaml")
# Estimate Genetic Map
map <- est_map(cross=magic, maxit=10, tol=1e-6, error_prob=0.02, map_function = "morgan", cores=2, quiet=F, save_rf=T)
magic$gmap <- map
max(map$`1`)  # for chrom 1
max(map$`2`)  # for chrom 2
max(map$`3`)  # for chrom 3
max(map$`4`)  # for chrom 4
max(map$`5`)  # for chrom 5
max(map$`6`)  # for chrom 6
max(map$`7`)  # for chrom 7
max(map$`8`)  # for chrom 8
max(map$`9`)  # for chrom 9
max(map$`10`)  # for chrom 10
magic$gmap <- map
saveRDS(map, file="example.Rdata")
linkage_map <- qtl2convert::map_list_to_df(map)
capture.output(map, file = "example.csv")
magic$gmap <- map

# Test scans
pr <- calc_genoprob(magic, map, error_prob=0.002, cores=4)
apr <- genoprob_to_alleleprob(pr)
kinship <- calc_kinship(pr)
grid <- calc_grid(magic$gmap, step=1)
pg <- scan1(pr, magic$pheno, kinship)
pg <- scan1(pr, magic$pheno, kinship, cores=4)
bin_pheno <- apply(magic$pheno, 2, function(a) as.numeric(a > median(a)))
find_peaks(pg, map, threshold=4.5, drop=1.5)


