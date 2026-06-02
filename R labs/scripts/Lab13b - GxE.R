###################################################
# CS/HS 541	- Plant breeding methods
# Lab 14 - GxE
# Roberto Fritsche-Neto
# roberto.neto@ncsu.edu
# Latest update: June 2, 2026
###################################################

#################################### MET analysis #########################
library(breedR)

#############################
# Load phenotypic data
#############################
pheno <- readRDS("pheno")
head(pheno)


###################################################
# Classical G×E using Mixed Models
###################################################
# Create genotype × environment identifier
# Each level represents one genotype in one environment
# first, we need to create a col for the GxN interaction
pheno$GN <- paste0(pheno$gid, pheno$N)
head(pheno)

###################################################
# Obtain trial design information
###################################################
# Number of replications
n.rep <- length(unique(pheno$rep))
# Number of environments
n.env <- length(unique(pheno$N))
n.rep
n.env

# Fitting genotype by environment models - with a common variance (diagonal model)
sol <- remlf90(fixed = SDM ~ N + rep, 
                random = ~ gid + GN,
                method = "em",
                data = pheno) 

# heritability and importance of GxE
sol$var
(h2g.plot <- sol$var[1,] / (sol$var[1,] + sol$var[2,] + sol$var[3,]))
(h2g <- sol$var[1,] / (sol$var[1,] + sol$var[2,]/loc + sol$var[3,]/(loc*blocks)))
(Hgxe <- sol$var[2,] / (sol$var[1,] + sol$var[2,] + sol$var[3,]))
sol$fit$AIC

# GxE deviations
head(sol$ranef$GN[[1]])


## Including genomics into the model - only the additive as an example
Ga <- readRDS("Ga")
Za <- model.matrix(~ -1 + gid, data = pheno)
colnames(Za) <- gsub("gid", "", colnames(Za), fixed = T)

sol2 <- remlf90(fixed = SDM ~ N + rep, 
                   random = ~ GN,
                   generic = list(GEBV = list(Za, Ga)),
                   method = "em",
                   data = pheno)

# heritability and importance of GxE
sol2$var
(h2g.plot.2 <- sol2$var[2,] / (sol2$var[1,] + sol2$var[2,] + sol2$var[3,]))
(h2g.2 <- sol2$var[2,] / (sol2$var[1,] + sol2$var[2,]/loc + sol2$var[3,]/(loc*blocks)))
(Hgxe.2 <- sol2$var[1,] / (sol2$var[1,] + sol2$var[2,] + sol2$var[3,]))
sol2$fit$AIC
# GxE deviations
head(sol2$ranef$GN[[1]])

## Including genomics and GxE into the model
Zge <- model.matrix(~ -1 + GN, data = pheno)
(N.levels <- diag(2))
colnames(N.levels) <- rownames(N.levels) <- unique(pheno$N)
N.levels
GxN <- kronecker(N.levels, Ga)
GxN[1:6, 1:6]
dim(GxN)

sol3 <- remlf90(fixed = SDM ~ N + rep, 
                   random = ~ 1,
                   generic = list(Ga = list(Za, Ga),
                                  GN = list(Zge, precision = GxN)),
                   method = "em",
                   data = pheno)

# heritability and importance of GxE
sol3$var
(h2g.plot.3 <- sol3$var[1,] / (sol3$var[1,] + sol3$var[2,] + sol3$var[3,]))
(h2g.3 <- sol3$var[1,] / (sol3$var[1,] + sol3$var[2,]/loc + sol3$var[3,]/(loc*blocks)))
(Hgxe.3 <- sol3$var[2,] / (sol3$var[1,] + sol3$var[2,] + sol3$var[3,]))
sol3$fit$AIC

# GxE deviations
head(sol3$ranef$GN[[1]])

# fitness of tested models
AIC.table <- data.frame(
  Model = c(
    "Classical GxE",
    "GBLUP",
    "Genomic GxE"
  ),
  AIC = c(
    sol$fit$AIC,
    sol2$fit$AIC,
    sol3$fit$AIC
  )
)

AIC.table # the smaller the better

###################################################
# Environmental covariance matrix
# Identity matrix assumes: environments are independent
###################################################
# If you want to include ERM and GRM into your model
#GxE <- kronecker(ERM, Ga)

################################################################################
############## Stability and adaptability - Finlay-Wilkinson ###################
################################################################################
library(statgenGxE)

## Create a TD object from dropsPheno.
dropsTD <- statgenSTA::createTD(data = pheno, 
                                genotype = "gid", 
                                trial = "N")

## Perform a Finlay-Wilkinson analysis for all trials.
dropsFW <- gxeFw(TD = dropsTD, 
                 trait = "SDM")
summary(dropsFW)

# let's take a look at the output
names(dropsFW)
dropsFW$estimates
dropsFW$envEffs
dropsFW$fittedGeno

## Create line plot for Finlay Wilkinson analysis.
plot(dropsFW, plotType = "line")

################################################################################
############################# GGE-Biplot Analysis ##############################
################################################################################
library(metan)

model.gge <- gge(pheno, N, gid, SDM, svp = "symmetrical")

model.gge$SDM

(a <- plot(model.gge, type = 1)) # basic plot
(b <- plot(model.gge, type = 2)) # Mean performance vs. stability
(c <- plot(model.gge, type = 3)) # Which-won-where
######## the end ##########