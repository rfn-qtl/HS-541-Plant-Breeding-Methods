#############################################
# CS/HS 541	- Plant breeding methods
# Lab 11 - GWAS
# Roberto Fritsche-Neto
# roberto.neto@ncsu.edu
# Latest update: June 1, 2026
#############################################
#-------------------------------------------------------------------------------
# Load Libraries and Setup Environment
#-------------------------------------------------------------------------------
# Loading all dependencies at the start ensures missing packages are caught early.
library(bigmemory)
library(biganalytics)
library(compiler)
library(PCAtools)
library(knitr)
library(breedR)
library(SNPRelate)
library(tidyverse) # For advanced data wrangling (tidyr) and plotting (ggplot2)
library(gaston)
library(CMplot)
library(reshape2)
library(ggplot2)

# Source external toolkits for GAPIT and FarmCPU engines
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")

#-------------------------------------------------------------------------------
# Data Ingestion and Quality Control
#-------------------------------------------------------------------------------
# Load phenotypes (2-step adjusted phenotypes, e.g., dBLUPs)
Y <- readRDS("pheno2step")
head(Y)
dim(Y)

# loading the marker file
M <- readRDS("M")
dim(M)

# Sync and filter genotypes: Filter and sort the marker matrix to perfectly match 
# the rows and order of the phenotype file's genetic IDs (gid)
M <- M[match(Y$gid, rownames(M)),]
dim(M)
head(M[,1:6])
tail(M[,1:6])
all(Y$gid == rownames(M))

# Load physical map metadata (HapMap)
hapmap <- readRDS("hapmap")
head(hapmap)
colnames(hapmap) <- c("SNP", "Chrm", "Pos")
dim(hapmap)
str(hapmap)

# Double check that the columns of matrix M exactly match the map rows
all(colnames(M) == hapmap$SNP)

# Create a structure matching GAPIT/FarmCPU requirements (Taxa column + markers)
M2 <- data.frame(Taxa = rownames(M), M)
head(M2[,1:6])
dim(M2)

#-------------------------------------------------------------------------------
# Population Structure Assessment via PCA
#-------------------------------------------------------------------------------
# We look for population stratification using the additive genomic relationship 
# matrix (Ga) to calculate principal components, which will serve as fixed covariates.
Ga <- readRDS("Ga")
# Reorder matrix rows and columns to match phenotype IDs
Ga <- Ga[match(Y$gid, rownames(Ga)),match(Y$gid, rownames(Ga))]
Ga[1:5, 1:5]
dim(Ga)

# Perform Principal Component Analysis (PCA)
metadata <- data.frame(row.names = colnames(Ga))
p <- pca(Ga, metadata = metadata, removeVar = 0.0)
# Structural diagnostic plots
screeplot(p, axisLabSize = 18, titleLabSize = 22)
biplot(p)
# Extract PC loading coordinates to pass as covariates (Covariate Variables)
PCA <- p$loadings
dim(PCA)
head(PCA[,1:3])

#-------------------------------------------------------------------------------
# 3. Formatting Datasets for FarmCPU
#-------------------------------------------------------------------------------
# myY : Column 1 = Taxa/GID, Column 2 = Phenotype (dBLUPs used as weights aren't built-in)
head(Y)
myY <- Y[,c(1, 3)] # using dBLUPs because there is no space to add weights in this package
head(myY)

myGD <- M2     # Genomic Data matrix
myGM <- hapmap # Genomic Map metadata
myCV <- PCA    # Covariate matrix (Population Structure)

#-------------------------------------------------------------------------------
# Multi-Locus GWAS: Running FarmCPU
#-------------------------------------------------------------------------------
# FarmCPU controls False Positives/Negatives by using a Fixed Effect Model 
# and a Random Effect Model iteratively. 

# Model A: Baseline model tracking kinship only
model.K <- FarmCPU(Y = myY, 
                   GD = myGD, 
                   GM = myGM, 
                   cutOff = 0.05, 
                   p.threshold = 0.05/nrow(myGM), 
                   MAF.calculate = T)
head(model.K$GWAS)
# To bypass the package in terms of kinship, you can create a diagonal matrix and feed the model

# Model B: Structural correction using 1 Principal Component
model.PC1 <- FarmCPU(Y = myY, 
                     GD = myGD, 
                     GM = myGM, 
                     CV = myCV[,1], 
                     cutOff = 0.05, 
                     p.threshold = 0.05/nrow(myGM), 
                     MAF.calculate = T)

# Model C: Structural correction using 2 Principal Components
model.PC2 <- FarmCPU(Y = myY, 
                     GD = myGD, 
                     GM = myGM, 
                     CV = myCV[,1:2], 
                     cutOff = 0.05, 
                     p.threshold = 0.05/nrow(myGM), 
                     MAF.calculate = T)

#-------------------------------------------------------------------------------
# Empirical Significance Thresholding (Permutations)
#-------------------------------------------------------------------------------
# Bonferroni corrections can be overly conservative due to Linkage Disequilibrium.
# We determine an empirically adjusted P-value threshold using 10 permutations.
# usually we use more than 100 permutations

p.adj <- FarmCPU.P.Threshold(Y = myY, 
                             GD = myGD, 
                             GM = myGM, 
                             trait = "dblups", 
                             theRep = 10)

p.Bonf <- 0.05 / ncol(M) # Standard Bonferroni comparison
print(paste("Empirical Adjusted P-threshold:", p.adj))
print(paste("Bonferroni P-threshold:", p.Bonf))

# Transforming thresholds to -log10 scale
(thr.p.adj <- -log10(p.adj))
(thr.p.Bonf <- -log10(p.Bonf))

# Re-run FarmCPU using our new optimized empirical threshold
model.p.adj <- FarmCPU(Y = myY, 
                       GD = myGD, 
                       GM = myGM, 
                       CV = myCV[,1:2], 
                       p.threshold = p.adj, 
                       threshold.output = p.adj, 
                       MAF.calculate = TRUE, 
                       cutOff = c((p.adj*ncol(M)), 0.05))

#-------------------------------------------------------------------------------
# Evaluation and Cross-Model Comparisons
#-------------------------------------------------------------------------------
# Merging GWAS P-value outputs from across our structural models
gwas.all <- data.frame(model.K$GWAS[,1:4], 
                       model.PC1$GWAS[, 4], 
                       model.PC2$GWAS[, 4], 
                       model.p.adj$GWAS[, 4])

colnames(gwas.all)[4:7] <- c("K", "1PC", "2PC", "p.adj")
head(gwas.all)

# Evaluate Pearson correlation matrix across models to evaluate control inflation
cor(gwas.all[,4:7], use = "pairwise.complete.obs")

# Generate Circular Manhattan and Q-Q Plot
CMplot(gwas.all, type = "p", plot.type = c("c", "q"),
       cir.chr.h=1.5,
       signal.line=1,
       width=10, height=10,
       chr.labels = paste("chrmo", " ", unique(hapmap$Chrm), sep=""),
       threshold=c(p.adj, 0.05/ncol(M)),threshold.col = c('black','red'),
       multracks = TRUE, chr.den.col = NULL, file.output = TRUE,
       signal.col = NULL, outward = TRUE, amplify = T, ylab = colnames(gwas.all[4:7]))

#-------------------------------------------------------------------------------
# Post-GWAS Analytics: QTL Effects & Single-Marker Metrics
#-------------------------------------------------------------------------------
# Extract significantly associated SNPs/QTLs
qtl <- model.p.adj$GWAS[which(model.p.adj$GWAS$P.value < p.adj), ] 
print(qtl)

# Visualization of Phenotype distribution relative to Minor Allele Copies (ggplot2)
if (nrow(qtl) > 0) {
  plot_data <- data.frame(dBLUP = Y$dBLUP)
  # Dynamic extraction of markers matching significant QTLs
  qtl_markers <- as.data.frame(M[, qtl$SNP, drop = FALSE])
  plot_data <- cbind(plot_data, qtl_markers)
  
  # Melt long formatting for clean multi-panel facet plotting
  plot_long <- plot_data %>% 
    pivot_longer(cols = -dBLUP, names_to = "SNP", values_to = "Copies") %>% 
    mutate(Copies = as.factor(round(Copies)))
  
  g_box <- ggplot(plot_long, aes(x = Copies, y = dBLUP, fill = Copies)) +
    geom_boxplot(outlier.alpha = 0.5, alpha = 0.7) +
    facet_wrap(~SNP, scales = "free_x") +
    labs(x = "Number of Minor Allele Copies", y = "SDM Phenotype (dBLUP)", 
         title = "Phenotypic Variance Across Genomic Classes") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2")
  print(g_box)
}


# Print summary tables of allele counts per class
for(i in seq_along(qtl$SNP)) {
  cat("\n--- SNP Allele Distribution Table:", qtl$SNP[i], "---\n")
  print(kable(table(round(M[, qtl$SNP[i]]))))
}

# Single-marker variance partitioning estimations
# Va = 2 * p * q * alpha^2
qtl$Va <- round(2 * qtl$maf * (1 - qtl$maf) * qtl$effect^2, 3)
qtl$ha <- round(qtl$Va / var(Y$dBLUP), 2)  # Proportion of phenotypic variance explained
qtl$Ac <- round(sqrt(qtl$ha), 2)           # Marker accuracy
qtl$R2 <- qtl$ha^2                         # Coefficient of determination proxy
print(kable(qtl, caption = "Identified QTL Genetic Architecture Attributes"))

#-------------------------------------------------------------------------------
# Joint-QTL Analysis via Mixed Model Framework
#-------------------------------------------------------------------------------
# If multiple QTLs exist, single-marker estimates may overlay or overestimate effects. 
# We build a comprehensive framework tracking population structure + joint mixed components.
fit.lm.data <- data.frame(dblups = Y$dBLUP,
                          gid = M2[,1],
                          PC1 = PCA[,1], 
                          PC2 = PCA[,2],
                          M[,as.character(qtl$SNP)])

# Normalize column naming conventions for modeling step
snp_cols_idx <- 5:ncol(fit.lm.data)
colnames(fit.lm.data)[snp_cols_idx] <- paste0("SNP", 1:nrow(qtl))

fit.lm.data$gid <- as.factor(fit.lm.data$gid)
fit.lm.data <- droplevels(fit.lm.data)

# Matrix mapping design framework
Za <- model.matrix(~ -1 + gid, data = fit.lm.data)

# Run Restricted Maximum Likelihood (REML) multi-locus marker regression
# Fixed effects handle structure (PCs) while targets are treated as random variances
fit.mlm <- remlf90(fixed = dblups ~ PC1 + PC2, # add the number of PCs you used in the GWAS
                   random = ~ SNP1, # add all the significant SNPs in the model
                   generic = list(Ga = list(Za, Ga)), 
                   method = "em",
                   data = fit.lm.data)

# Compute Joint-Heritability components
var.comp <- fit.mlm$var
ha.snp <- var.comp[1] / sum(var.comp)
print(paste("Joint Narrow-sense Heritability explained by SNP1:", round(ha.snp, 4)))

# Predict Genomic Estimated Breeding Values (GEBVs) exclusively based on validated QTLs
bv <- as.matrix(M[, qtl$SNP, drop = FALSE]) %*% qtl$effect
qtl_accuracy <- cor(bv, Y$dBLUP)
print(paste("Predictive accuracy of QTL-based breeding values:", round(qtl_accuracy, 4)))

#-------------------------------------------------------------------------------
# Linkage Disequilibrium (LD) Decay and Local Architecture Mapping
#-------------------------------------------------------------------------------
# Chrm by Chrm
head(hapmap)
hapmap <- hapmap[hapmap$Chrm == qtl$Chrm,]
M <- M[, hapmap$SNP]
dim(M); dim(hapmap)

map1 <- data.frame(snp.id = as.integer(1:dim(hapmap)[1]), 
                   snp.rs.id = hapmap$SNP, 
                   chr = as.integer(hapmap$Chrm), 
                   pos = as.integer(hapmap$Pos))
head(map1)

# Build a standardized Genomic Data Structure (GDS) file payload
snpgdsCreateGeno(gds.fn = "./toy.gds",                              # gds filename
                 genmat = M,                                        # markers matrix
                 sample.id = rownames(M),                           # individual names
                 snp.id = map1$snp.id,                               # snp id (deve ser único)
                 snp.rs.id = map1$snp.rs.id,                         # snp name (pode não ser único)
                 snp.chromosome = map1$chr,                          # cromossomo
                 snp.position = map1$pos,                            # posição dentro do cromossomo (bp)
                 #snp.allele = map$allele,                           # alelos (ref / alt)
                 snpfirstdim = FALSE)                                # argumento para matriz n x s

# Loading gds
genofile <- snpgdsOpen(filename = "./toy.gds")
genofile

# Estimate pairwise Linkage Disequilibrium ($r^2$) matrix
LDs <- abs(as.matrix(snpgdsLDMat(genofile, method = "corr", slide = 0, num.thread = 8)$LD))
dim(LDs)
rownames(LDs) <- colnames(LDs) <- hapmap$SNP
LDs[1:5, 1:5]

# reshape the LD file to a long format
lds <- melt(LDs)
head(lds)
dim(lds)

# function to estimate the distance between two markers A nd B
naive_dist <- function(A, B){
  result = matrix(ncol = length(B), nrow = length(A))
  for (i in 1:length(A))
    for (j in 1:length(B))
      result[i,j] = abs(A[i] - B[j])
  result
}

# estimating the distance between markers
distance <- naive_dist(as.matrix(hapmap$Pos), as.matrix(hapmap$Pos)) 
rownames(distance) <- colnames(distance) <- hapmap$SNP
distance2 <-  melt(distance)
head(distance2)
dim(distance2)

# let's remove some heavy files
rm(distance)

# and combine the information
pairLD <- data.frame(marker1 = lds$Var1, 
                     marker2 = lds$Var2, 
                     r2 = lds$value, 
                     dist = distance2$value)
head(pairLD)

# let's remove some heavy files
rm(distance2); rm(lds);  

# Plot fine-mapped local Linkage Disequilibrium surrounding the focal QTL
qtl_idx <- which(hapmap$SNP == qtl$SNP[1])
window_range <- max(1, qtl_idx - 5):min(nrow(hapmap), qtl_idx + 5)

LD.plot(LDs[window_range, window_range], 
        hapmap$Pos[window_range],
        graphical.par = list(cex = 1.3, bg = "gray"), 
        polygon.par = list(border = NA))


# LD Decay plot ans estimates
# Create the plot with a smoothed trendline (LOESS or Generalized Additive Model)
# If the dataset is massive (millions of rows), sub-sample for plotting
set.seed(123)
ld_sub <- pairLD[sample(1:nrow(pairLD), 10000),]

# Strip diagonal comparisons
ld_sub <- ld_sub[ld_sub$dist != 0,] 
dim(ld_sub)

# estimate the average LD
(averageLD <- round(mean(na.omit(ld_sub$r2)), 2))

# finally, the famous plot
p_draft <- ggplot(ld_sub, aes(x = dist, y = r2)) +
  geom_point(alpha = 0.1, color = "gray") +  # Light points to show density
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "red", size = 1.2) + # Decay line
  labs(
    title = "Linkage Disequilibrium (LD) Decay",
    x = "Physical Distance (kb)",
    y = expression(LD~(r^2))
  ) +
  theme_minimal() +
  ylim(0, 1)

# Extract the calculated smooth line data
plot_data <- ggplot_build(p_draft)$data[[1]]

# Find the x (distance) value where y (R2) is closest to 0.2
# plot_data$x maps to distance_kb, plot_data$y maps to R2
distance_intercept <- approx(x = plot_data$y, y = plot_data$x, xout = 0.2)$y

# Print the result
print(paste("LD decays to r2 = 0.2 at:", round(distance_intercept, 2), "kb"))

# finally, the famous plot
ggplot(ld_sub, aes(x = dist, y = r2)) +
  geom_point(alpha = 0.1, color = "gray") +  # Light points to show density
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "red", size = 1.2) + # Decay line
  labs(
    title = "Linkage Disequilibrium (LD) Decay",
    x = "Physical Distance (kb)",
    y = expression(LD~(r^2))
  ) +
  
  # Add the horizontal line (e.g., at r2 = 0.2)
  geom_hline(yintercept = averageLD, linetype = "dashed", color = "blue", size = 0.8) +
  
  # Add the vertical line (e.g., at 150 kb)
  geom_vline(xintercept = distance_intercept, linetype = "dashed", color = "blue", size = 0.8) +
  
  theme_minimal() +
  ylim(0, 1)

################### the end ############