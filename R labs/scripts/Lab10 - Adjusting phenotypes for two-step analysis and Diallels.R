#############################################
# CS/HS 541	- Plant breeding methods
# Lab 10 - Adjusting phenotypes and Diallels
# Roberto Fritsche-Neto
# roberto.neto@ncsu.edu
# Latest update: July 8, 2025
#############################################

##################### adjusting phenotypes for two-set analysis ###############
# loading files
pheno <- readRDS("pheno")
head(pheno)

# creating a col for interaction or nested effects
pheno$GN <- as.factor(paste0(pheno$gid, pheno$N))
head(pheno)

# adjusting the model 
library(breedR)

# First way - genotypes as fixed
solF <- remlf90(fixed = SDMadj ~ N + gid, 
                random = ~ rep + GN,
                method = "em",
                data = pheno) 
                
(sol.fixed <- solF$fixed$gid)

# GE deviations
head(solF$ranef$GN[[1]])

# Second way - deregressing phenotypes
solR <- remlf90(fixed = SDMadj ~ N, 
                random = ~ rep + gid + GN, 
                method = "em", 
                data = pheno)
(sol.random <- solR$ranef$gid)

# realiability
(rel <- 1 - (sol.random[[1]][,2])^2/solR$var[2])
mean(rel)

#accuracy
(Ac.MPS <- sqrt(mean(rel)))

# derregressing BLUPS
dblup <- sol.random[[1]][,1]/rel
dblup

# weights
c <- 0.5 # this constant is amount of variation that we expect not be explained by markers
(hg <- round(solR$var[2] / (solR$var[2] + solR$var[3]),2)) # heritability
(w <- (1 - hg)/(c + (1 - rel) / rel*hg))

# combining the data
pheno2step <- data.frame(gid = rownames(sol.fixed[[1]]), 
                      BLUE = sol.fixed[[1]][,1], 
                      dBLUP = dblup, 
                      w = w)
head(pheno2step)

# comparing the adjustments
par(mfrow = c(1,2))
hist(pheno2step$BLUE, col = "red", main = "BLUE", xlab = NULL)
hist(pheno2step$dBLUP, col = "blue", main = "dBLUP", xlab = NULL)
dev.off()

cor(pheno2step$BLUE, pheno2step$dBLUP)
plot(pheno2step$BLUE, pheno2step$dBLUP, xlab = "BLUES", ylab = "dBLUP" )
abline(lm(pheno2step$dBLUP ~ pheno2step$BLUE), col = "red")

# selecting only the hybrids
(hybrids <- unique(pheno[pheno$type == "sc","gid"]))
pheno2step <- droplevels.data.frame(pheno2step[pheno2step$gid %in% hybrids,])

# saving the newest file
saveRDS(pheno2step, "pheno2step")

################################ MATING DESIGNS #########################

# loading data and a library
pheno <- readRDS("pheno")
A <- readRDS("A")
Ga <- readRDS("Ga")
Gd <- readRDS("Gd")

# creating some incidence matrices for genetic effects
# additive
Za <- model.matrix(~ -1 + gid, data = pheno)
Za[1:5,1:5]
dim(Za)
colnames(Za) <- gsub("gid", "", colnames(Za), fixed = T)
all(pheno$gid %in% colnames(Za))
# dominance
pheno$gidD <- pheno$gid
Zd <- model.matrix(~ -1 + gidD, data = pheno)
Zd[1:5,1:5]
dim(Zd)
colnames(Zd) <- colnames(Za)
all(pheno$gid %in% colnames(Zd)) 

#################### parents and hybrids together ####################
library(breedR)

## model I (parents and hybrids non-related)
modI <- remlf90(SDMadj ~ N + rep,
                 random = ~gid,
                 data = pheno)

(varI <- modI$var)
(hg.I <- varI[1] / (varI[1] + varI[2]))


## model A (pedigree)
modA <- remlf90(fixed = SDMadj ~ N + rep, 
                random = ~ 1, 
                generic = list(A = list(Za, A)), 
                method = "em",
                data = pheno)
(varA <- modA$var)
(hg.A <- varA[1] / (varA[1] + varA[2]))

## model Ga
modGa <- remlf90(fixed = SDMadj ~ N + rep, 
                      random = ~ 1, 
                      generic = list(Ga = list(Za, Ga)), 
                      method = "em",
                      data = pheno)

(varGa <- modGa$var)
(hg.Ga <- varGa[1] / (varGa[1] + varGa[2]))

## model Ga & Gd
modGaGd <- remlf90(fixed = SDMadj ~ N + rep, 
                       random = ~ 1,
                       generic = list(Ga = list(Za, Ga),
                                      Gd = list(Zd, precision = Gd)),
                       method = "em",
                       data = pheno)

(varGaGd <- modGaGd$var)
(hg.GaGd <- varGaGd[1] / (varGaGd[1] + varGaGd[2] + varGaGd[3]))


###################
# BLUPS
##################
# joint the data
BLUPS <- data.frame(gid = colnames(Za),
                    I = modI$ranef$gid[[1]][,1],
                    A = modA$ranef$A[[1]][,1],
                    Ga = modGa$ranef$Ga[[1]][,1],
                    GaGd = modGaGd$ranef$Ga[[1]][,1])

head(BLUPS)

require(PerformanceAnalytics)
# correlarion betwwen BLUPS ontaines from differente methods
chart.Correlation(BLUPS[,-1], histogram = TRUE, pch = 1, method = "pearson")


############### Half-diallel or NCII - GS-based method ########################
# selecting only the data that correspond to hybrids
phenoSC <- droplevels.data.frame(pheno[pheno$type == "sc",])

# Ga for females
Ga[1:7, 1:7]
Ga.f <- Ga[1:7, 1:7]
# Ga for males
Ga[8:14, 8:14]
Ga.m <- Ga[8:14, 8:14]
# crosses
Gd[1:4, 1:4]
Gd.sc <- Gd[15:ncol(Gd), 15:ncol(Gd)]

# creating incidence matrices for females and males
head(pheno)
Za.f <- model.matrix(~ -1 + female, data = phenoSC)
Za.f[1:5,1:5]
dim(Za.f)
colnames(Za.f) <- gsub("female", "", colnames(Za.f), fixed = T)
all(phenoSC$female %in% colnames(Za.f))

Za.m <- model.matrix(~ -1 + male, data = phenoSC)
Za.m[1:5,1:5]
dim(Za.m)
colnames(Za.m) <- gsub("male", "", colnames(Za.m), fixed = T)
all(phenoSC$male %in% colnames(Za.m))

Zd.sc <- model.matrix(~ -1 + gid, data = phenoSC)
Zd.sc[1:5,1:5]
dim(Zd.sc)
colnames(Zd.sc) <- gsub("gid", "", colnames(Zd.sc), fixed = T)
all(phenoSC$gid %in% colnames(Zd.sc))

## model
NCII <- remlf90(fixed = SDMadj ~ N + rep, 
                   random = ~ 1,
                   generic = list(Ga.f = list(Za.f, Ga.f),
                                  Ga.m = list(Za.m, precision = Ga.m),
                                  Gd.sc = list(Zd.sc, precision = Gd.sc)),
                   method = "em",
                   data = phenoSC)

# components of variance
(varNCII <- NCII$var)

# heritabilities
(ha.f <- varNCII[1] / (varNCII[1] + varNCII[2] + varNCII[3] + varNCII[4]))
(ha.m <- varNCII[2] / (varNCII[1] + varNCII[2] + varNCII[3] + varNCII[4]))
(hd <- varNCII[3] / (varNCII[1] + varNCII[2] + varNCII[3] + varNCII[4]))
(ha <- (varNCII[1] + varNCII[2])/ 
    (varNCII[1] + varNCII[2] + varNCII[3] + varNCII[4]))
(hg <- (varNCII[1] + + varNCII[2] + varNCII[3])/ 
    (varNCII[1] + varNCII[2] + varNCII[3] + varNCII[4]))

# breeding values
BV.f <- NCII$ranef$Ga.f[[1]]
BV.m <- NCII$ranef$Ga.m[[1]]
head(BV.f)
head(BV.m) 

# CGA
GCA.f <- BV.f / 2
GCA.m <- BV.m / 2
head(GCA.f) 
head(GCA.m)   

#SCA 
SCA <- NCII$ranef$Gd.sc[[1]]
head(SCA)

################################################################################
###### Building heterotic pools (dataset does not help, but let's anyway) ######
################################################################################
pheno <- readRDS("pheno")
# create a col for the interaction female x males
pheno$fm <- paste0(pheno$female, pheno$male) 
head(pheno)

# For that we need to run a full diallel using G matrices
# Let's create the SCA matrix - only the parents
Ga <- readRDS("Ga")
Gp <- Ga[1:14, 1:14]
Gd.mf <- kronecker(Gp, Gp)
dim(Gd.mf)
Gd.mf[1:5, 1:5] # SCA
# and incidence matrix for them
Zf <- model.matrix(~ -1 + female, data = pheno)
Zm <- model.matrix(~ -1 + male, data = pheno)
Zp <- Zf + Zm
Zp[1:10, 1:4]
colnames(Zp) <- gsub("female", "", colnames(Zp), fixed = T)
dim(Zp) 
Zp[1:10, 1:4]

# colnames for Gd.mf
aux <- expand.grid(colnames(Gp), colnames(Gp))
head(aux)
aux$fm <- paste0(aux$Var1, aux$Var2)
colnames(aux) <- c("female", "male", "fm") 
head(aux)
dim(aux)
colnames(Gd.mf) <- rownames(Gd.mf) <- aux$fm
Gd.mf[1:5, 1:5]
all(pheno$fm %in% colnames(Gd.mf))

# now, th incidence matrix
Z.sca.seen <- model.matrix(~ -1 + fm, data = pheno)
Z.sca.seen[1:5, 1:5]
dim(Z.sca.seen)
colnames(Z.sca.seen) <- gsub("fm", "", colnames(Z.sca.seen), fixed = T)

Z.sca.missed <- matrix(0, 
                       ncol = (ncol(Gd.mf) - ncol(Z.sca.seen)), 
                       nrow = nrow(Z.sca.seen))
Z.sca.missed[1:5, 1:5]
dim(Z.sca.missed)
colnames(Z.sca.missed) <- colnames(Gd.mf)[-match(colnames(Z.sca.seen), colnames(Gd.mf))]
setdiff(aux$fm, unique(pheno$fm))

# combine both matrices
Z.sca <- cbind(Z.sca.seen, Z.sca.missed)
dim(Z.sca)

# finally, run the model
FD <- remlf90(fixed = SDMadj ~ N + rep, 
                random = ~ 1,
                generic = list(GCA = list(Zp, Gp),
                               SCA = list(Z.sca, precision = Gd.mf)),
                method = "em",
                data = pheno)

FD$ranef$GCA
FD$ranef$SCA

# merge datasets
SCA <- data.frame(fm = colnames(Z.sca), 
                  SCA = FD$ranef$SCA[[1]]$value)
SCA <- merge(SCA, aux)
head(SCA)
tail(SCA)
dim(SCA)
SCA <- SCA[,2:4]

library(reshape2)
d <- reshape(SCA, idvar = "female", timevar = "male", direction = "wide")
dim(d)
d <- d[,-1]
colnames(d) <- gsub("SCA.", "", colnames(d), fixed = T)
d <- d[, match(colnames(d), sort(as.numeric(colnames(d))) )]
rownames(d) <- colnames(d)
diag(d) <- 0
d[1:5, 1:5]

# clustering in heterotic pools
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

k <- kmeans(scale(d), centers = 2, nstart = 25)
fviz_cluster(k, data = scale(d), main = "Heterotic Pools")
k$size
sum(k$size)
(heterotic.pools <- data.frame(sort(k$cluster)))

# Selecting testers
GCA <- FD$ranef$GCA[[1]]
rownames(GCA) <- colnames(Gp)
GCA

testers <- data.frame()
for(i in 1:2){
  HP <- GCA[rownames(heterotic.pools)[heterotic.pools$sort.k.cluster. == i],]
  parent <- rownames(HP)[which.max(HP$value)]
  testers <- rbind(testers,
                   data.frame(
                     HP = i,
                     Tester = parent
                   ))}

testers
########## the end ####################