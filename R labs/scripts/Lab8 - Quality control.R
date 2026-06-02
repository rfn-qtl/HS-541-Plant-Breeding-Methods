#######################################################
# CS/HS 541	- Plant breeding methods
# Lab 8 - Quality control of marker and phenotypic data
# Roberto Fritsche-Neto
# roberto.neto@ncsu.edu
# Latest update: Jun 1, 2026
#######################################################

############################################################
######################### Phenotypes #######################
############################################################

# loading the phenotyping file
pheno <- read.table("../data/pheno.txt", header = T, na.strings = "NA")
head(pheno) # show the first six rows
tail(pheno) # show the last six rows
str(pheno) # present the flies' structure

# adjusting to factors
pheno[,1:9] <- lapply(pheno[,1:9], factor)
str(pheno)

# X and Y coordinate as a variate
pheno$X <- as.numeric(pheno$row)    
pheno$Y <- as.numeric(pheno$col)
str(pheno)

# phenotypic correlation between traits
cor_matrix <- round(cor(pheno[,10:12], use = "pairwise.complete.obs"), 2)
cor_matrix

library("PerformanceAnalytics")
chart.Correlation(
  as.matrix(na.omit(pheno[,10:12])),
  histogram = TRUE,
  pch = 16
)

############################################################
# Outlier detection
############################################################
# We fit a linear model accounting for design effects and
# identify observations with unusually large residuals.
boxplot(pheno$SDM, col = "red")
#install.packages("lme4")
library(lme4)
# outlier detection and elimination
fit <- lm(SDM ~ type + row + col + N + gid, data = pheno)
library(car)
(outlier <- names(outlierTest(fit)$p))
pheno[outlier,  "SDM"]
pheno[outlier, "SDM"] <- NA


# another way, based on Cook’s Distance Diagnostics and Outlier Detection 
#The common cutoff for defining outliers is 4/(n−p), 
# where n is the number of observations and p is the number of parameters
library(remotes)
library(datacooks)
df1 <- datacooks(fit, threshold= 4, clean= FALSE)
dim(df1)
df2= datacooks(fit, threshold= 4, clean= TRUE)
dim(df2)

############################################################
# Experimental design and spatial distribution
############################################################
library(desplot)
d1 <- desplot(pheno, SDM ~ X*Y, out1 = rep, out2 = N,
              out2.gpar=list(col = "green", lwd = 1, lty = 1))
print(d1)

############################################################
# Normality assessment
############################################################
# Many statistical models assume normally distributed
# residuals rather than normally distributed raw data.
# Nevertheless, evaluating the trait distribution provides
# useful information about potential transformations.
# First lets check using patterns
shapiro.test(rnorm(length(pheno$SDM))) # normal distribution
shapiro.test(runif(length(pheno$SDM))) # uniform distribution
# then, 
shapiro.test(pheno$SDM)

# Apply automatic normalization methods and select the
# transformation that best approximates normality.
#install.packages("bestNormalize")
require(bestNormalize)
SDMadj <- bestNormalize(pheno$SDM, standardize = FALSE, allow_orderNorm = TRUE, out_of_sample = FALSE)
SDMadj$chosen_transform
shapiro.test(SDMadj$x.t)
pheno$SDMadj <- SDMadj$x.t
head(pheno)

############################################################
# Residual diagnostics
############################################################
# Quartile‐Quartile (Q‐Q) normality plot for residuals
fit <- lm(SDM ~ type + row + col + N + gid, data = pheno)
fit2 <- lm(SDMadj ~ type + row + col + N + gid, data = pheno)

par(mfrow = c(2,2)) # organize the plot window in 1 row and 2 col
qqnorm(resid(fit))
qqline(resid(fit), col = "red")
qqnorm(resid(fit2))
qqline(resid(fit2), col = "blue")
hist(pheno$SDM, col = "red", main = "SDM", xlab = "SDM")
hist(pheno$SDMadj, col = "blue", main = "Adjusted SDM", xlab = "Adjusted SDM")
dev.off()

# saving the newest pheno file
str(pheno)
head(pheno)
saveRDS(pheno, "pheno")


###############################################################################
################################## Markers ###################################
###############################################################################
# Marker QC removes:
# 1. SNPs with excessive missing data
# 2. Rare alleles (low MAF)
# 3. Poorly genotyped markers
# 4. Markers requiring imputation
# This improves genomic relationship estimates and
# downstream prediction accuracy.

#loading file
geno <- readRDS("../data/geno")
head(geno)
dim(geno)

# let's check the SNPs and genotyping errors
(reads <- unique(unlist(apply(geno[,5:19], 2, unique))))

# create the hapmap file
hapmap <- geno[,c(1, 3, 4)]
head(hapmap)
tail(hapmap)
dim(hapmap)
str(hapmap)
hapmap$chrom <- as.factor(hapmap$chrom)
hapmap$pos <- as.numeric(hapmap$pos)
str(hapmap)

# marker file
M <- t(geno[,6:19])
M[1:5, 1:5]
dim(M)

# check if they are equal
all(colnames(M) == rownames(hapmap)) 

##################################### Quality control ######################
#BiocManager::install("impute")
#devtools::install_github(repo = 'italo-granato/snpReady', ref = 'dev')

library(ASRgenomics)
#aux2 <- ASRgenomics::snp.recode(as.matrix(M))$Mrecode
#aux2[1:4, 1:5]

library(snpReady)
QC <- raw.data(data = as.matrix(M), 
               frame = "wide", 
               hapmap = hapmap, 
               sweep.sample = 1, #0.80
               maf = 0.15, # 0.05 # Minor allele frequency threshold
               # Rare variants contribute little information and may
               # inflate sampling variance
               call.rate = 0.90,   # Minimum sample call rate
               # Marker must be scored in at least 90% of individuals
               base = TRUE, 
               imput = TRUE, 
               imput.type = "wright", 
               outfile = "012", 
               plot = F)

# The report of the quality control approaches
QC$report

# get the newest dataset of markers
M <- QC$M.clean
dim(M)
M[1:5, 1:5]

# and the hapmap
hapmap <- QC$Hapmap
dim(hapmap)
head(hapmap)

############################################################
# Remove heterozygous markers
############################################################
# The population consists of highly inbred lines.
# Therefore, heterozygous calls may indicate:
# - genotyping errors
# - residual heterogeneity
# - sample contamination
# We retain only markers containing homozygous genotypes.
non.het.markers <- apply(M, 2, function(x){all(x != 1)})
sum(non.het.markers)
M <- M[, non.het.markers]
dim(M)

# then, correcting the map object again
hapmap <- hapmap[non.het.markers, ]
head(hapmap)
dim(hapmap)

# verifying if all the markers into M matrix are in the same order into map
identical(as.character(hapmap$rs), colnames(M))

############################################################
# Linkage disequilibrium (LD) pruning
############################################################
# Adjacent markers often carry redundant information
# because they are highly correlated.
# LD pruning removes markers that provide nearly identical
# genetic information, reducing:
# - computational burden
# - matrix collinearity
# - storage requirements
# Here, markers with r² > 0.98 are considered redundant.

# SNPRelate package
# to install the package
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("SNPRelate")
library(SNPRelate)

# hapmap
map2 <- data.frame(snp.id = as.integer(1:dim(hapmap)[1]), 
                   snp.rs.id = hapmap$rs, 
                   chr = as.integer(hapmap$chrom), 
                   pos = as.integer(hapmap$pos))
head(map2)
tail(map2)

# creating the GDS file
snpgdsCreateGeno(gds.fn = "./toy.gds",                              # gds filename
                 genmat = M,                                        # markers matrix
                 sample.id = rownames(M),                           # individual names
                 snp.id = map2$snp.id,                               # snp id 
                 snp.rs.id = map2$snp.rs.id,                         # snp name 
                 snp.chromosome = map2$chr,                          # cromossomo
                 snp.position = map2$pos,                            # positionon chrom (bp)
                 #snp.allele = map$allele,                           # alleles (ref / alt)
                 snpfirstdim = FALSE)                                # argument for matrix n x s


# Loading gds
genofile <- snpgdsOpen(filename = "./toy.gds")
genofile


# prune markers by MAF 0.10, CR 0.90 e LD 0.99 (r2)
snps_pruned <- snpgdsLDpruning(gdsobj = genofile,         # gds file
                               remove.monosnp = TRUE,     # only bialellic
                               maf = 0.15,                # MAF
                               missing.rate = 0.90,       # CR
                               method = "corr",           # r2 (method of LD)
                               slide.max.bp = 100000,     # slide window 
                               ld.threshold = 0.98,       # LD threshold
                               start.pos = "first")       # start at

# get SNP ids
snps_pruned <- unlist(unname(snps_pruned))
snps_pruned
length(snps_pruned) # good quality snps 

# removing SNPs prunned by LD from M matrix
M <- M[, as.numeric(snps_pruned)]
M[1:5, 1:5]
dim(M)

# correcting the map object
head(hapmap)
hapmap <- hapmap[as.numeric(snps_pruned), ] # removing SNPs pruned by LD
head(hapmap)
dim(hapmap)

# verifying if all the markers into M matrix are in the same order into map
identical(as.character(hapmap$rs), colnames(M))

# close GDS
snpgdsClose(genofile)

############################################################
# Create hybrid marker matrix
############################################################
# Single-cross hybrids are generated by averaging the
# parental marker scores.
# Example:
# Female = 0
# Male   = 2
# Hybrid = (0 + 2)/2 = 1
# which represents the expected heterozygous genotype.

phenoSC <- pheno[pheno$type == "sc",]
dim(phenoSC)
head(phenoSC)
sc.grid <- expand.grid(sort(unique(phenoSC$female)), sort(unique(phenoSC$male)))

M.female <- M[match(sc.grid[, 1], rownames(M)),]
M.male <- M[match(sc.grid[, 2], rownames(M)),]
dim(M.female); dim(M.male) ; dim(M)
M.SC <- (M.female + M.male) / 2
dim(M.SC)
rownames(M.SC) <- apply(sc.grid, 1, function(x){paste0(x[1], x[2])})
M.SC[1:5, 1:5]
M[1:5, 1:5]

# Finally combine them, parents and hybrids
M <- rbind(M, M.SC)
dim(M)
all(unique(phenoSC$gid) %in% rownames(M))

# saving files
saveRDS(hapmap, "hapmap")
saveRDS(M, "M")
################# the end #############