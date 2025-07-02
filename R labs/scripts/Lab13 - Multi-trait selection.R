###################################################
# CS/HS 541	- plant breeding methods
# Lab 13 - Multi-trait selection
# Roberto Fritsche-Neto
# roberto.neto@ncsu.edu
# Latest update: July 1, 2025
###################################################
# loading packages
require(foreach)
require(doParallel)
require(doMC)
library(car)
library(ggplot2)
library(breedR)

# setting the number of cores that will be used
detectCores()
registerDoParallel(cores = detectCores()-1) # type the number of cores you want to use
getDoParWorkers()

#loading phenotypes
pheno <- readRDS("pheno")
head(pheno)
dim(pheno)
str(pheno)

# reorganizing the file
pheno <- pheno[,c(1:9, 13:14, 10:12)]
head(pheno)

# index for traits
traits <- colnames(pheno)[12:14]
  
# estimate de phenotypic correlation
pheno.cor <-round(cor(pheno[,traits], use = "pairwise.complete.obs", method = "pearson"), 2)
corrplot::corrplot(pheno.cor, method = 'number', type = "lower", diag = F, col = c("red", "orange", "green", "blue"))

# reorganize de data
pheno.melted <- reshape2::melt(pheno, measure.vars = traits)
head(pheno.melted)  

# load the GRM
Ga <- readRDS("Ga")

# running all GS single-traits in parallel 
results.st <- foreach(i = 1:length(traits), 
                        .packages = c("breedR", "car"), 
                        .combine = "rbind",
                        .export = c("remlf90", "outlierTest"),
                        .multicombine = TRUE, 
                        .errorhandling = "remove",
                        .verbose = TRUE    
  ) %do% {

    # subset the data  
    smpl <- droplevels.data.frame(pheno.melted[pheno.melted$variable == traits[i],])
    
    # outlier detection and elimination
    fit <- lm(value ~ rep + gid + N + gid:N, data = smpl)
    outlier <- names(outlierTest(fit)$p)
    smpl[outlier, "value"] <- NA
    
    Za <- model.matrix(~ -1 + gid, data = smpl)
    colnames(Za) <- colnames(Ga)

    # MME using only the classical experimental design, and gid as random
    sol <- remlf90(fixed = value ~ 1, 
                   generic = list(BV = list(Za, Ga)), 
                   data = smpl, 
                   method = "em")
    
    # reorganizing BLUPS
    BLUPS <- data.frame(gid = colnames(Ga),
                        trait = traits[i],
                        GEBV = sol$ranef$BV[[1]][,1])

}

head(results.st)
tail(results.st)

# GEBVs per trait
GEBVs <- tapply(results.st$GEBV, list(results.st$gid, results.st$trait), mean)
head(GEBVs)

######################## single step (ss) multi-trait (MT) GS ################

# covariance matrices 
covg <- matrix(c(1, cov(pheno[,12], pheno[,13], use = "pairwise.complete.obs"), 
                 cov(pheno[,12], pheno[,13], use = "pairwise.complete.obs"), 1), 2, 2)
colnames(covg) <- rownames(covg) <- colnames(pheno)[12:13]
covg

covr <- matrix(c(1, 0, 0, 1), 2, 2)
colnames(covr) <- rownames(covr) <- colnames(pheno)[2:3]
covr

sol <- remlf90(fixed = cbind(pheno[,2], pheno[,3]) ~ N + rep, 
               generic = list(GEBV = list(Za, Ga, var.ini = covg)), 
               data = pheno, 
               method = "em", 
               var.ini = list(covr))

## predicted values for the test set
blupsA <- as.matrix(sol$ranef[[1]][[1]])
rownames(blupsA) <- colnames(Ga)
blupsB <- as.matrix(sol$ranef[[1]][[2]])
rownames(blupsB) <- colnames(Ga)

# heritabilities
(hm.A <- sol$var[[1]][1,1]/sum(sol$var[[1]][1,1], sol$var[[2]][1,1]))
(hm.B <- sol$var[[1]][2,2]/sum(sol$var[[1]][2,2], sol$var[[2]][2,2]))
(LL <- logLik(sol)[1])

GEBV.MT = data.frame(
  gid = colnames(Ga),
  GEBV_A = blupsA[,1], 
  EBV_B = blupsB[,1])

# G and R varcomp via MT-GBLUP
sol$var

# the new correlation between the traits
cor(GEBV.MT[,2:3])
# correlation among traits via MT-GBLUP
cov2cor(sol$var$GEBV)
# the residual correlation
cov2cor(sol$var$Residual)

# So, the only thing now is just weight the GEBV and obtain the SI
# define the economic weights per trait
ecoW <- c(1, 1) # in this case two times more for grain yield
# then, the vector of SI per genotype
MT <- as.matrix(GEBV.MT[,2:3]) %*% ecoW  
head(MT)

########################### Selection indices ##################################
# phenotypic covariance between traits
P <- cov(pheno[, 12:13], use = "pairwise", method = "pearson")

# genetic covariance between traits
G <- cov(GEBVs[,2:3])

# Smith-Hazel
# define the economic weights per trait
ecoW <- c(1, 1) # in this case two times more for grain yield
# then, the selection weights per trait
(b.sh <- solve(as.matrix(P)) %*% as.matrix(G) %*% as.matrix(ecoW))
# Finally, the vector of SI per genotype
SH <- GEBVs[,2:3] %*% b.sh  
head(SH)

# Pasek-Baker
# define the desired genetic gains per traits in genetic standard deviations
desired <- c(1, 1)
# G correlation matrix x desired genetic gains in standard deviations
(b.pb <- solve(cov(scale(GEBVs[,2:3]))) %*% as.matrix(desired))
PB <- as.matrix(scale(GEBVs[,2:3])) %*% b.pb
head(PB)

# correlation between GEBVs and selection indices
GEBVs <- as.data.frame(GEBVs[,2:3])
GEBVs$SH <- SH
GEBVs$PB <- PB
GEBVs$MT <- MT
head(GEBVs)
cor(GEBVs)

# let's see which method provides the best KPIs
# the biggest average and the smallest variation
apply(cor(GEBVs[, 3:5]), 2, mean) / apply(cor(GEBVs[, 3:5]), 2, sd)

# defining a threshold to select the individuals to advance
is <- 15 #15%
quantile(GEBVs$SH, ((100 - is) / 100))
GEBVs$Seleted <- GEBVs$SH > quantile(GEBVs$SH, ((100 - is) / 100))
# the number of pre-selected parents
sum(GEBVs$Seleted)
head(GEBVs)
# saving the file
write.csv(GEBVs, "GEBVS.csv")

# 3D plot showing the selected materials and parents
library(ggplot2)
library(ggpubr)
a <- ggplot(GEBVs, aes(x = SDM, 
                       y = SRA,
                       color = Seleted)) + geom_point()
a
#################### the end ####################################