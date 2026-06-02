###################################################
# CS/HS 541	- Plant breeding methods
# Lab 9 - Pedigree and Genomic Relationship Matrix
# Roberto Fritsche-Neto
# roberto.neto@ncsu.edu
# Latest update: Jun 1, 2026
###################################################

############################## PEDIGREE ###############################

# install.packages("pedigreemm")
library(pedigreemm) #load pedigreemm package

###################################################
# Load pedigree information
###################################################
# The pedigree file contains:
# gid    = genotype identifier
# male   = male parent
# female = female parent
# Pedigree information allows estimation of
# expected genetic relationships among individuals.
ped  <- read.table("../data/ped.txt", header = TRUE)
head(ped, 15)
tail(ped)
str(ped)

###################################################
# Pedigree ordering
###################################################
# Relationship matrices require ancestors to appear
# before descendants.
# editPed() automatically reorders individuals from
# oldest generations to newest generations.
args(editPed)
ped2 <- editPed(ped$male, ped$female, ped$gid)
head(ped2) # the last col is generation
tail(ped2)

###################################################
# Create pedigree object
###################################################
# The pedigree object stores all parent-offspring
# connections and serves as the basis for calculating
# additive genetic relationships.
args(pedigree)
ped3 <- pedigree(ped2$sire, ped2$dam, ped2$label)
head(ped3)
tail(ped3)

###################################################
# Numerator relationship matrix (A)
###################################################
# A represents expected additive genetic relationships.
# Examples:
# Relationship          Expected value
# ------------------------------------
# Individual with self       1.0
# Parent - offspring         0.5
# Full siblings              0.5
# Half siblings              0.25
#
# Diagonal elements:
#   self-relationship
#
# Off-diagonal elements:
#   pairwise additive relationship

A <- getA(ped3)
A <- as.matrix(A) 
dim(A)
A[1:14, 1:14] # the first seven parents
A[47:55, 47:55] # the last hybrids

# let`s save or A matrix`
saveRDS(A, "A")

###################################################
# Population structure analysis
###################################################
# PCA summarizes relationship patterns among individuals.
# Individuals with similar ancestry tend to cluster
# together in principal component space.
#
# This approach is widely used to:
# - identify subpopulations
# - detect family groups
# - identify outliers
# - evaluate breeding pools

# svd decomposition - by individuals
# Singular value decomposition (SVD) is mathematically
# equivalent to PCA and allows extraction of major
# axes of genetic variation.
svdG <- svd(A, nu = ncol(A), nv = nrow(A))
plot(cumsum((svdG$d[1:ncol(A)])^2/sum(svdG$d^2)), ylab = "proportion accumulated", xlab = "number of individuals", col = "red")

# obtainig the eigenvectors and eigenvalues
pcsG <- A %*% svdG$v
rownames(pcsG) <- colnames(pcsG) <- rownames(A)
dim(pcsG)
pcsG[1:14,1:5]

###################################################
# Variance explained
###################################################

var_explained <- round(
  svdG$d^2 / sum(svdG$d^2) * 100,
  2
)
head(var_explained)
cumsum(var_explained)

# 3D graph 
library(scatterplot3d)
scatterplot3d(pcsG[,1], pcsG[,2], pcsG[,3], xlab = "PC1", ylab = "PC2", zlab = "PC2", axis = TRUE, color = "red", highlight.3d = FALSE, box = TRUE, angle = 50)

# 2D graph
par(mfrow = c(1,2))
plot(x = pcsG[,1], y = pcsG[,2], xlab = "PC1", ylab = "PC2", col = "red", main = "PC 1 vs PC 2")
plot(x = pcsG[,1], y = pcsG[,3], xlab = "PC1", ylab = "PC2", col = "blue", main = "PC 1 vs PC 3")
dev.off()

###################################################
# Relationship heatmap
###################################################
# Heatmaps provide a visual representation of
# pairwise genetic relationships.
#
# Darker regions indicate closely related groups.
#
# Clusters often correspond to:
# - families
# - heterotic groups
# - breeding populations
#install.packages("superheat")

library(superheat)

superheat(A, pretty.order.rows = T, pretty.order.cols = T, col.dendrogram = T, clustering.method = "kmeans", 
          dist.method = "euclidean",  bottom.label.text.size = 2, left.label.text.size = 2, legend.text.size = 5)

# saving the graph for papers
png("heatmapA.png", width = 8, height = 8, res = 400, units = "in")
superheat(A, pretty.order.rows = T, pretty.order.cols = T, col.dendrogram = T, clustering.method = "kmeans", 
          dist.method = "euclidean",  bottom.label.text.size = 2, left.label.text.size = 2, legend.text.size = 5)
dev.off()

###################################################
# Pedigree-based inbreeding
###################################################
# Inbreeding measures the probability that two
# alleles are identical by descent.
#
# Higher values indicate increased homozygosity
# resulting from mating related individuals.
(Inbreeding <- inbreeding(ped3))

###################################################
# Genomic relationship matrix (GRM)
###################################################
# Unlike pedigree relationships, genomic
# relationships are estimated directly from DNA
# markers.
#
# Advantages:
# - captures Mendelian sampling
# - accounts for realized relationships
# - detects hidden relatedness
# - generally improves prediction accuracy

# read marker data
M  <- readRDS("M")
dim(M)
head(M[,1:6])
tail(M[,1:6])

#source("https://bioconductor.org/biocLite.R")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install("impute")
#devtools::install_github(repo = 'italo-granato/snpReady', ref = 'dev')
library(snpReady)
library(ASRgenomics)
# creates the K matrix - genomic relatioship matrix 
args(G.matrix)

G <- snpReady::G.matrix(M, method = "VanRaden", format = "wide", plot = F)

############################### additive #######################
Ga <- G$Ga
dim(Ga)

# let`s load our A matrix and check if they are the same population
A <- readRDS("A")
G2A <- match.G2A(A = A, G = Ga, clean = TRUE, ord = F, mism = TRUE,
                 RMdiff = TRUE)
Ga <- G2A$Gclean
dim(Ga)
Ga[1:4, 1:4]
Ga[52:55, 52:55]

###################################################
# Matrix rescaling
###################################################

# Pedigree and genomic matrices often operate on
# different numerical scales.
# Rescaling facilitates visual comparison between
# A and G but does not alter relationship rankings.

library(scales)
rescale.G <- function(x){
  out <- scales::rescale(c(x), to = c(0, 2)) 
  out <- matrix(out, nrow = nrow(x), byrow = T)
  rownames(out) <- colnames(out) <- colnames(x)
  return(out)
}

Ga <- rescale.G(Ga)
dim(Ga)
range(Ga)
A[1:4, 1:4]
Ga[1:4, 1:4]
A[52:55, 52:55]
Ga[52:55, 52:55]

###################################################
# Relationship diagnostics - Ga
###################################################
# Diagnostics identify:
# - unusually high self-relatedness
# - duplicated samples
# - possible labeling errors
# - contaminated genotypes
check_Ga <- kinship.diagnostics(K = Ga, diagonal.thr.small = 0.8,
                                diagonal.thr.large = 1.2, duplicate.thr = 0.95)
check_Ga$plot.diag
check_Ga$plot.offdiag

###################################################
# Dominance relationship matrix
###################################################
# Ga:
# additive genetic relationships
# Gd:
# dominance genetic relationships
#
# Additive effects:
#   transmitted to offspring
# Dominance effects:
#   interaction between alleles at the same locus
#
# Dominance is particularly important in hybrid
# breeding because it contributes to heterosis.
Gd <- G$Gd[colnames(Ga), colnames(Ga)]
dim(Gd)
Gd[1:4, 1:4]
Gd[52:55, 52:55]

rescale.Gd <- function(x){
  out <- scales::rescale(c(x), to = c(0, 1)) 
  out <- matrix(out, nrow = nrow(x), byrow = T)
  rownames(out) <- colnames(out) <- colnames(x)
  return(out)
}

Gd <- rescale.Gd(Gd)
dim(Gd)
Gd[1:4, 1:4]
Gd[52:55, 52:55]
range(Gd)

###################################################
# Relationship diagnostics - Gd
###################################################

check_Gd <- kinship.diagnostics(K = Gd, diagonal.thr.small = 0.8,
                                diagonal.thr.large = 1.2, duplicate.thr = 0.95)
check_Gd$plot.diag
check_Gd$plot.offdiag

#################### # graphs analysis ##############################

# svd decomposition for Ga and Gd
Ga_pca <- kinship.pca(K = Ga, ncp = 14, label = T, ellipses = T)
Ga_pca$plot.pca
Ga_pca$plot.scree

Gd_pca <- kinship.pca(K = Gd, ncp = 14, label = T, ellipses = T)
Gd_pca$plot.pca
Gd_pca$plot.scree

###################################################
# Population structure from markers
###################################################
# PCA based on genomic relationships often reveals
# finer population structure than pedigree records.
#
# It can identify:
# - hidden ancestry
# - subpopulations
# - misclassified genotypes
# - heterotic groups

kinship.heatmap(K = Ga, dendrogram = TRUE, row.label = T,
                col.label = T)

kinship.heatmap(K = Gd, dendrogram = TRUE, row.label = T,
                col.label = T)

###################################################
# Inbreeding and effective population size
###################################################
# Genomic inbreeding:
# F = diagonal(G) - 1
# Positive values indicate excess homozygosity.
#
# Effective population size (Ne) quantifies the
# amount of genetic diversity maintained within
# the breeding population.
#
# Small Ne:
#   - increased drift
#   - increased inbreeding
#   - reduced genetic gain potential
#
# Large Ne:
#   - greater long-term diversity
#   - better sustainability of selection
# Inbreeding per individual - considering just the parents

(Fi <- round(diag(Ga)-1, 2)[1:14])
# Effective size per individual
(Ne.i <- 1/(2*Fi))
# Ne of population
(Ne.pop <- sum(Ne.i))
# Population endogamy rate
(Fi.pop <- 1/(2*Ne.pop))

# for a more accurate estimate of F and Ne check 
# Caballero et al. (2022)
# https://doi.org/10.1186/s12711-022-00772-0

###################################################
# Pedigree versus genomic relationships
###################################################
corA_G <- cor(c(A), c(Ga), use = "complete.obs")
cat("\nCorrelation between A and G:", round(corA_G, 3), "\n")
# A high correlation indicates agreement between
# pedigree expectations and marker-derived
# relationships.
#
# Differences arise because genomic relationships
# capture realized Mendelian segregation whereas
# pedigree relationships represent only expected
# relationships.

############################################################
# Single-step relationship matrix (H)
############################################################
# The H matrix combines:
# A = pedigree-based relationships
# G = marker-based relationships
#
# Advantages:
# - Uses all available information
# - Incorporates genotyped and non-genotyped individuals
# - More accurate than pedigree alone
# - Basis of modern single-step genomic prediction
# References:
# Aguilar et al. (2010)
# Legarra et al. (2009)

############################################################
# Construct H
############################################################
library(ASRgenomics)

H <- H.matrix(
  A = A,
  G = solve(Ga),
)

str(H)
dim(H)
H[1:5,1:5]

#saving the final kinship versions
saveRDS(Ga, "Ga")
saveRDS(Gd, "Gd")
saveRDS(H,  "H")
############# the end ##############