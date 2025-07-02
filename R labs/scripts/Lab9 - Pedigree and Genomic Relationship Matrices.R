###################################################
# CS/HS 541	- Plant breeding methods
# Lab 9 - Pedigree and Genomic Relationship Matrix
# Roberto Fritsche-Neto
# roberto.neto@ncsu.edu
# Latest update: Jun 30, 2025
###################################################

############################## PEDIGREE ###############################

# install.packages("pedigreemm")
library(pedigreemm) #load pedigreemm package

# read in the pedigree data
ped  <- read.table("../data/ped.txt", header = TRUE)
head(ped, 15)
tail(ped)
str(ped)

#editPed function orders the pedigre from oldest to newest	
args(editPed)
ped2 <- editPed(ped$male, ped$female, ped$gid)
head(ped2) # the last col is generation
tail(ped2)

#constructs the pedigree object
args(pedigree)
ped3 <- pedigree(ped2$sire, ped2$dam, ped2$label)
head(ped3)
tail(ped3)

# creates the A matrix (relationship matrix, which means 2 times de kinship matrix)
A <- getA(ped3)
A <- as.matrix(A) 
dim(A)
A[1:14, 1:14] # the first seven parents
A[47:55, 47:55] # the last hybrids

# let`s save or A matrix`
saveRDS(A, "A")

####################
# graphs analysis
####################

# svd decomposition - by individuals
svdG <- svd(A, nu = ncol(A), nv = nrow(A))
plot(cumsum((svdG$d[1:ncol(A)])^2/sum(svdG$d^2)), ylab = "proportion accumulated", xlab = "number of individuals", col = "red")

# obtainig the eigenvectors and eigenvalues
pcsG <- A %*% svdG$v
rownames(pcsG) <- colnames(pcsG) <- rownames(A)
dim(pcsG)
pcsG[1:14,1:5]

# PCA graphs
# proportion explained by the first componentes
axispcs <- paste((round(svdG$d[1:ncol(A)]^2/sum(svdG$d^2)*100))[1:3], "%", sep = "")
axispcs

# 3D graph 
library(scatterplot3d)
scatterplot3d(pcsG[,1], pcsG[,2], pcsG[,3], xlab = axispcs[1], ylab = axispcs[2], zlab = axispcs[3], axis = TRUE, color = "red", highlight.3d = FALSE, box = TRUE, angle = 50)

# 2D graph
par(mfrow = c(1,2))
plot(x = pcsG[,1], y = pcsG[,2], xlab = axispcs[1], ylab = axispcs[2], col = "red", main = "PC 1 vs PC 2")
plot(x = pcsG[,1], y = pcsG[,3], xlab = axispcs[1], ylab = axispcs[3], col = "blue", main = "PC 1 vs PC 3")
dev.off()

##############
# heatmaps
##############
#install.packages("superheat")
library(superheat)

superheat(A, pretty.order.rows = T, pretty.order.cols = T, col.dendrogram = T, clustering.method = "kmeans", 
          dist.method = "euclidean",  bottom.label.text.size = 2, left.label.text.size = 2, legend.text.size = 5)

# saving the graph for papers
png("heatmapA.png", width = 8, height = 8, res = 400, units = "in")
superheat(A, pretty.order.rows = T, pretty.order.cols = T, col.dendrogram = T, clustering.method = "kmeans", 
          dist.method = "euclidean",  bottom.label.text.size = 2, left.label.text.size = 2, legend.text.size = 5)
dev.off()

(Inbreeding <- inbreeding(ped3))


############################## GRM ###############################
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

# rescale the matrix to match with the expected values
library(scales)
rescale.G <- function(x){
  out <- scales::rescale(c(x), to = c(0, 2)) 
  out <- matrix(out, nrow = nrow(x), byrow = T)
  rownames(out) <- colnames(out) <- colnames(x)
  return(out)
}

Ga <- rescale.G(Ga)
dim(Ga)
Ga[1:4, 1:4]
Ga[52:55, 52:55]

# let's check the Ga matrix
check_Ga <- kinship.diagnostics(K = Ga, diagonal.thr.small = 0.8,
                                diagonal.thr.large = 1.2, duplicate.thr = 0.95)
check_Ga$plot.diag
check_Ga$plot.offdiag

########################### dominance kinship #######################
Gd <- G$Gd[colnames(Ga), colnames(Ga)]
dim(Gd)
Gd[1:4, 1:4]
Ga[52:55, 52:55]


rescale.Gd <- function(x){
  out <- scales::rescale(c(x), to = c(0, 1)) 
  out <- matrix(out, nrow = nrow(x), byrow = T)
  rownames(out) <- colnames(out) <- colnames(x)
  return(out)
}

Gd <- rescale.Gd(Gd)
dim(Gd)
Gd[1:4, 1:4]
Ga[52:55, 52:55]

# let's check the Gd matrix
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

############### heatmaps ##############

kinship.heatmap(K = Ga, dendrogram = TRUE, row.label = T,
                col.label = T)

kinship.heatmap(K = Gd, dendrogram = TRUE, row.label = T,
                col.label = T)

#saving the final kinship versions
saveRDS(Ga, "Ga")
saveRDS(Gd, "Gd")

############################ Inbreeding and Ne ########################
# Inbreeding per individual - considering just the parents
(Fi <- round(diag(Ga)-1, 2)[1:14])
# Effective size per individual
(Ne.i <- 1/(2*Fi))
# Ne of population
(Ne.pop <- sum(Ne.i))
# Population endogamy rate
(Fi.pop <- 1/(2*Ne.pop))

############# the end ##############