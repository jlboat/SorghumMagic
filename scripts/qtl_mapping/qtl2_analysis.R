library(qtl2)

magic <- ("~/Desktop/MAGIC1")
# Read cross using newly estimated map
magic <- read_cross2("~/Desktop/MAGIC1/magic.yaml")
map <- insert_pseudomarkers(magic$gmap, step=1)
pr <- calc_genoprob(magic, map, error_prob=0.002, cores=4)
apr <- genoprob_to_alleleprob(pr)

# Estimate Kinship
kinship <- calc_kinship(pr)
grid <- calc_grid(magic$gmap, step=1)
pr_grid <- probs_to_grid(pr, grid)
kinship_grid <- calc_kinship(pr_grid)

# Scan
pg <- scan1(pr, magic$pheno, kinship, cores=4)
bin_pheno <- apply(magic$pheno, 2, function(a) as.numeric(a > median(a)))
find_peaks(pg, map, threshold=4.5, drop=1.5)
plot(pg,map, lodcolumn=1, col= "slateblue")
plot(pg,map, lodcolumn=2, col= "violetred", add = TRUE)
plot(pg,map, lodcolumn=3, col= "red", add = TRUE)
plot(pg,map, lodcolumn=4, col= "black", add = TRUE)
# Permutation Scan @5000 permutations
operm <- scan1perm(pr, magi$pheno, n_perm=5000)
summary(operm)
##For adding legends to the QTL graphs
legend(3,22, legend=c("PCR1CU","PCR2CU", "PCR1FL", "PCR2FL"), col=c("slateblue","violetred", "red", "black"),lty=1, cex=0.8)



