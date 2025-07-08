###################################################
# CS/HS 541	- Plant breeding methods
# Lab 12 - Genomic Selection
# Roberto Fritsche-Neto
# roberto.neto@ncsu.edu
# Latest update: July 8, 2025
###################################################

# loading data
Y <- readRDS("pheno2step")
head(Y)
str(Y)
dim(Y)

# add a column to model the dominance
Y$gidD <- Y$gid
head(Y)

# loading kernels
Ga <- readRDS("Ga")
dim(Ga)
Gd <- readRDS("Gd")
dim(Gd)

# selecting just the hybrids
Ga <- Ga[match(Y$gid, colnames(Ga)), match(Y$gid, colnames(Ga))]
dim(Ga)
all(Y$gid == colnames(Ga))
Gd <- Gd[match(Y$gid, colnames(Gd)), match(Y$gid, colnames(Gd))]
dim(Gd)
all(Y$gid == colnames(Gd))

# creating the Za for additive effects
Za <- model.matrix(~ -1 + gid, data = Y)
colnames(Za) <- gsub("gid", "", colnames(Za), fixed = T)
dim(Za)
Za[1:6, 1:6]

# creating the Zd for additive effects
Zd <- model.matrix(~ -1 + gidD, data = Y)
colnames(Zd) <- gsub("gidD", "", colnames(Za), fixed = T)
dim(Zd)
Zd[1:6, 1:6]

##################### GS using A model ####################
library(breedR)

# Cross-validation settings 
nfold <- 5 # training and validation sizes
replicates <- 5 # repeat this procedure 5 times

output_A <- data.frame() # create an empty file to storage KPI
GEBVs <- data.frame() # create an empty file to GEBV
pred.vs.obs <- data.frame()

for (r in 1:replicates) {
  cat("Processing the replicate", r, "\n")
  ## test sets
  sets <- sample(rep(1:nfold, length.out = n <- dim(Y)[1])) # different orders of the same vector
  
  for (k in 1:nfold) {
    cat("Processing the fold", k, "\n")
    
    # creating the training and testing sets
    itest = which(sets == k)
    itrain = which(sets != k)
    YGS <- Y
    YGS[itest, "dBLUP"] <- NA # it creates the validation set
    
    # Running the model
    sol <- remlf90(fixed = dBLUP ~ 1, 
                   generic = list(A = list(Za, Ga)), 
                   data = YGS, 
                   method = "em")
  
    # reorganizing BLUPS
    BLUPS <- data.frame(gid = colnames(Za),
                        GEBV = sol$ranef$A[[1]][,1])
    predicted.blups <- BLUPS[BLUPS$gid %in% YGS$gid[itest], ]
    observed.blues <- Y[Y$gid %in% YGS$gid[itest],]
    observed.blues <- observed.blues[match(predicted.blups$gid, observed.blues$gid),]
      
    output_A <- rbind(output_A, data.frame(
      method = "A",
      rep = r,
      fold = k,
      h2m = round(sol$var[1,]/sum(sol$var), 2),
      PA = round(cor(predicted.blups[, 2], observed.blues[, 3]), 2),
      AIC = sol$fit$AIC,
      LRT = sol$fit$`-2logL`
      )
    )
    
    pred.vs.obs <- rbind(pred.vs.obs, data.frame(
      predicted.blups,
      observed.blues
      )
    )
    
    # storage genetic breeding estimated values
    GEBVs <- rbind(GEBVs, BLUPS)
    
  }
}


# estimate the KPI mean
head(output_A)
tail(output_A)
(resultA <- apply(output_A[,4:6], 2, mean)) 

# and looking at the distribution
library(ggplot2)
ggplot(output_A, aes(x = method, y = PA)) + 
  geom_boxplot(notch = TRUE, outlier.colour="red", outlier.shape=16, outlier.size=2) +
  geom_jitter(position=position_jitter(0.2))

# lets' see the overall relationship between predicted and observed values
head(pred.vs.obs)
dim(pred.vs.obs)

# obtain the average
pred.vs.obsA <- data.frame(
  gid = names(tapply(pred.vs.obs$GEBV, pred.vs.obs$gid, mean)),
  BV = tapply(pred.vs.obs$GEBV, pred.vs.obs$gid, mean),
  BLUES = tapply(pred.vs.obs$dBLUP, pred.vs.obs$gid, mean)
)

# the accuracy and distribution
library(plotly)
cor(pred.vs.obsA[,2], pred.vs.obsA[,3])
(p <- plot_ly(pred.vs.obsA, 
              x = ~BV, 
              y = ~BLUES, 
              type = 'scatter', 
              mode = 'markers',  
              text = ~paste('gid: ', gid)))

# and the GEBV to compare RR-BLUP and GBLUP
dim(GEBVs)
head(GEBVs)

# lets organize the breeding values and pick the mean 
GEBVs <- as.matrix(tapply(GEBVs$GEBV, GEBVs$gid, mean))
head(GEBVs)
GEBVs <- GEBVs[(match(colnames(Ga), rownames(GEBVs))),]
GEBVs <- as.matrix(GEBVs)
colnames(GEBVs) <- "GEBV"
head(GEBVs)
dim(GEBVs)

##################### GS using A + D model + BLUES and weights ####################

# Cross-validation settings 
nfold <- 5 # training and validation sizes
replicates <- 5 # repeat this procedure 5 times

output_AD <- data.frame() # create an empty file to storage KPI
GEGVs <- data.frame() # create an empty file to GEBV

for (r in 1:replicates) {
  cat("Processing the replicate", r, "\n")
  ## test sets
  sets <- sample(rep(1:nfold, length.out = n <- dim(Y)[1])) # different orders of the same vector
  
  for (k in 1:nfold) {
    cat("Processing the fold", k, "\n")
    
    # creating the training and testing sets
    itest = which(sets == k)
    itrain = which(sets != k)
    YGS <- Y
    YGS[itest, "BLUE"] <- NA # it creates the validation set
    
    # Running the model
    sol <- remlf90(fixed = BLUE ~ 1, 
                   generic = list(A = list(Za, Ga),
                                  D = list(Zd, precision = Gd)),
                   weights = YGS$w,
                   data = YGS, 
                   method = "em")
    
    # reorganizing BLUPS
    BLUPS <- data.frame(gid = colnames(Za),
                        GEBV = sol$ranef$A[[1]][,1],
                        GEGV = sol$ranef$A[[1]][,1] + sol$ranef$D[[1]][,1]
                        )
    predicted.blups <- BLUPS[BLUPS$gid %in% YGS$gid[itest], ]
    observed.blues <- Y[Y$gid %in% YGS$gid[itest],]
    observed.blues <- observed.blues[match(predicted.blups$gid, observed.blues$gid),]
    
    output_AD <- rbind(output_AD, data.frame(
      method = "AD",
      rep = r,
      fold = k,
      h2m = round(sol$var[1,]/sum(sol$var), 2),
      PA = round(cor(predicted.blups[, 3], observed.blues[, 2]), 2),
      AIC = sol$fit$AIC,
      LRT = sol$fit$`-2logL`)
      )
    
    # retrieving only the genetic deviations
    blups <- data.frame(gid = colnames(Za),
                        A = sol$ranef$A[[1]][,1],
                        D = sol$ranef$D[[1]][,1],
                        GV = sol$ranef$A[[1]][,1] + sol$ranef$D[[1]][,1])
    
    GEGVs <- rbind(GEGVs, blups)
    
  }
}

# estimate the KPI mean
head(output_AD)
(resultAD <- apply(na.omit(output_AD[,4:7]), 2, mean)) 

# and looking at the PA distribution
library(ggplot2)
ggplot(output_AD, aes(x = method, y = PA)) + 
  geom_boxplot(notch = TRUE, outlier.colour = "red", outlier.shape = 16, outlier.size = 2)+
  geom_jitter(position=position_jitter(0.2))

########### the end ################