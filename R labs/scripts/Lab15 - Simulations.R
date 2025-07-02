###################################################
# CS/HS 541	- Plant breeding methods
# Lab 15 - Simulations
# Roberto Fritsche-Neto
# roberto.neto@ncsu.edu
# Latest update: July 2, 2025
###################################################
library(AlphaSimR)

#################### Parameters ##############
nQTL <- 350
add <- 0.22  
chip.size.1 <- 910
segSites <- chip.size.1 + nQTL
replicates <- 5
r <- 1:replicates
cycles <- 3

# Heritabilities by breeding stage
#        F2,   F3,   F4,   F5,   F6,   F7
h2 <- c(0.03, 0.15, 0.40, 0.60, 0.70, 0.80)
H2 <- c(0.06, 0.20, 0.45, 0.63, 0.72, 0.81)

# crossing block
n.parents <- 40
n.crosses <- 40
progenie.size <- 100

###########################
# crop history of evolution
###########################
set.seed(123)

history <- runMacs(nInd = 400, 
                   nChr = 7,
                   segSites = segSites / 7, 
                   inbred = TRUE,
                   species = "GENERIC", 
                   split = NULL, 
                   ploidy = 2L)

# then, add the simulation parameters to the simulated population
SP = SimParam$new(history)

# Lets add a trait, considering gamma distribution for QTL effect;  A + D effects
SP$restrSegSites(nQTL/7, chip.size.1/7)
SP$addTraitAD(nQtlPerChr = nQTL / 7, gamma = TRUE, mean = 0, var = 1, meanDD = add, varDD = 0.5)

# add a SNP chip to our population
SP$addSnpChip(chip.size.1/7)

# after the historical population (mimic the crop evolution, it is time to simulate the first years of a breeding program, like as a warm up process)
# setting the heritability of the trait, narrow and broad sense
SP$setVarE(h2 = h2[5], H2 = H2[5])
firstPop <- newPop(rawPop = history, isDH = TRUE, simParam = SP)

#########################################################
# simulating the founders
#########################################################
# select the first parents, full inbreed
set.seed(123)

founders <- selectInd(firstPop, nInd = n.parents, trait = 1, use = "pheno", gender = "B",
                      selectTop = TRUE, returnPop = TRUE, candidates = NULL,
                      simParam = SP)

#############################################################
# simulating 3 generations of trad breeding == burn-in
#############################################################
cat("Burn-in", "\n")

# intensity of selection
(nSelF3 <- n.crosses*progenie.size*0.1)
(nSelF4 <- nSelF3*0.1)
(nSelF5 <- nSelF4*0.1)
(nSelF6 <- 1)

pop.trad <- founders
perform <- data.frame()

for (i in 1:3) {
  cat("Processing the cycle", i, "\n")
  
  # Breeding Crosses (Year 1)
  f1 <- randCross(pop = pop.trad, nCrosses = n.crosses, nProgeny = 1, simParam = SP)
  
  # Single Plant Selection (Year 2)
  SP$setVarE(h2 = h2[1], H2 = H2[1])
  f2 <- self(pop = f1, nProgeny = progenie.size, simParam = SP)
  f2sel <- selectInd(pop = f2, nInd = nSelF3, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
  
  # Progeny Row Stage 1 (Year 3)
  SP$setVarE(h2 = h2[2], H2 = H2[2])
  f3 <- self(pop = f2, nProgeny = 1, simParam = SP)
  f3sel <- selectInd(pop = f3, nInd = nSelF4, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
  
  # Preliminary Yield Test (Year 4)
  SP$setVarE(h2 = h2[3], H2 = H2[3])
  f4 <- self(pop = f3sel, nProgeny = 1, simParam = SP)
  f4sel <- selectInd(pop = f4, nInd = nSelF5, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
  
  # Advanced Test (URN - Year 5)
  SP$setVarE(h2 = h2[4], H2 = H2[4])
  f5 <- self(pop = f4sel, nProgeny = 1, simParam = SP)
  f5sel <- selectInd(pop = f5, nInd = nSelF6, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
  
  # Advanced Yield (CA - Years 6-7)
  SP$setVarE(h2 = h2[5], H2 = H2[5])
  f6 <- self(pop = f5sel, nProgeny = 1, simParam = SP)
  f6sel <- selectInd(pop = f6, nInd = 1, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
  
  # Foundation Increase (Year 8)
  SP$setVarE(h2 = h2[6], H2 = H2[6])
  variety <- self(pop = f6sel, nProgeny = 1, simParam = SP)
  
  # then, select the best one regardless the family
  newparents <- f4
  pop.trad <- newparents
  
  # output
  perform <- rbind(perform, data.frame(population = as.numeric(meanG(pop.trad)), 
                                       cultivar = as.numeric(max(variety@gv))))
}

perform

#########################
# creating the first and small TS to be updated
# a population of 1152 - 3 chips
#########################
cat("TS", "\n")
#cross everyone
SP$setVarE(h2 = h2[4], H2 = H2[4])
TS0 <- randCross(pop = pop.trad, nCrosses = 24, nProgeny = 48, simParam = SP)

# recombining the offspring
TS00 <- newPop(rawPop = TS0, simParam = SP)

# obtaining full inbred lines by DH tech
TS000 <- makeDH(pop = TS00, nDH = 1, useFemale = TRUE, simParam = SP)
markers <- RRBLUP2(pop = TS000, traits = 1, use = "pheno", snpChip = 1, simParam = SP)

#Evaluate accuracy via an independent pop
Testing_set <- randCross(pop = pop.trad, nCrosses = 24, nProgeny = 48, simParam = SP)
Testing_set = setEBV(Testing_set, markers, simParam=SP)
Testing_set@ebv
(Ac.GS <- as.numeric(cor(gv(Testing_set), ebv(Testing_set))))
capture.output(c(Ac.GS = Ac.GS), file = "GS.accuracies.txt")

# storage the first results
# results C0
resultsC0 <-  data.frame(
  method = c("Trad", "Trad+GS"),
  rep = 1,
  cycle = 0,
  PM = as.numeric(meanG(pop.trad)),
  Va = as.numeric(varA(pop.trad)),
  Ac = c(as.numeric(sqrt(varG(pop.trad)/varP(pop.trad))), Ac.GS),
  Variety = perform[3,2],
  Years = c(5)
)


##############################################################
# simulating more cycles of breeding by the Traditional Method
##############################################################
cat("Trad", "\n")

out.Trad <- data.frame()

for (k in 1:length(r)) {
  cat("Processing Replicate:", k, "\n")
  
  pop.base <- pop.trad
  
  for (i in 1:cycles) {
    cat("Processing cycle", i, "Trad", "\n")
    # Breeding Crosses (Year 1)
    f1 <- randCross(pop = pop.base, nCrosses = n.crosses, nProgeny = 1, simParam = SP)
    
    # Single Plant Selection (Year 2)
    SP$setVarE(h2 = h2[1], H2 = H2[1])
    f2 <- self(pop = f1, nProgeny = progenie.size, simParam = SP)
    f2sel <- selectInd(pop = f2, nInd = nSelF3, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    # Progeny Row Stage 1 (Year 3)
    SP$setVarE(h2 = h2[2], H2 = H2[2])
    f3 <- self(pop = f2, nProgeny = 1, simParam = SP)
    f3sel <- selectInd(pop = f3, nInd = nSelF4, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    # Preliminary Yield Test (Year 4)
    SP$setVarE(h2 = h2[3], H2 = H2[3])
    f4 <- self(pop = f3sel, nProgeny = 1, simParam = SP)
    f4sel <- selectInd(pop = f4, nInd = nSelF5, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    # Advanced Test (URN - Year 5)
    SP$setVarE(h2 = h2[4], H2 = H2[4])
    f5 <- self(pop = f4sel, nProgeny = 1, simParam = SP)
    f5sel <- selectInd(pop = f5, nInd = nSelF6, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    # Advanced Yield (CA - Years 6-7)
    SP$setVarE(h2 = h2[5], H2 = H2[5])
    f6 <- self(pop = f5sel, nProgeny = 1, simParam = SP)
    f6sel <- selectInd(pop = f6, nInd = 1, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    # Foundation Increase (Year 8)
    SP$setVarE(h2 = h2[6], H2 = H2[6])
    variety <- self(pop = f6sel, nProgeny = 1, simParam = SP)
    
    # updating the parents and pop
    newparents <- f4
    pop.base <- newparents
    
    # output
    out.Trad <- rbind(out.Trad, data.frame(
      method = "Trad",
      rep = k,
      cycle = i,
      PM = as.numeric(meanG(pop.base)),
      Va = as.numeric(varA(pop.base)),
      Ac = as.numeric(sqrt(varG(pop.base)/varP(pop.base))),
      Variety = as.numeric(gv(variety)),
      Years = 5
    ))
  }
}

###################################################################
# simulating more cycles of breeding by the Traditional Method + GS
###################################################################
cat("Trad+GS", "\n")

out.Trad_GS <- data.frame()
pop.base <- pop.trad

for (k in 1:length(r)) {
  cat("Processing Replicate:", k, "\n")
  
  pop.base <- pop.trad
  TS <- TS000
  popList <- list(TS)
  snps <- markers
  
  for (i in 1:cycles) {
    cat("Processing cycle", i, "Trad+GS", "\n")
    
    f1 <- randCross(pop = pop.base, nCrosses = n.crosses, nProgeny = 1, simParam = SP)
    
    SP$setVarE(h2 = h2[1], H2 = H2[1])
    f2 <- self(pop = f1, nProgeny = progenie.size, simParam = SP)
    f2sel <- selectInd(pop = f2, nInd = nSelF3, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    SP$setVarE(h2 = h2[2], H2 = H2[2])
    f3 <- self(pop = f2sel, nProgeny = 1, simParam = SP)
    
    # GS
    f3.gs <- setEBV(pop = f3, solution = snps, value = "bv", simParam = SP)
    f3sel <- selectInd(pop = f3.gs, nInd = nSelF4, trait = 1, "ebv", gender = "B", selectTop = TRUE, simParam = SP)
    Ac.GS <- as.numeric(cor(gv(f3.gs), ebv(f3.gs)))
    
    # setting the heritability of the trait, narrow and broad sense
    SP$setVarE(h2 = h2[3], H2 = H2[3])
    f4 <- self(pop = f3sel, nProgeny = 1, simParam = SP)
    f4sel <- selectInd(pop = f4, nInd = nSelF5, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    # setting the heritability of the trait, narrow and broad sense
    SP$setVarE(h2 = h2[4], H2 = H2[4])
    f5 <- self(pop = f4sel, nProgeny = 1, simParam = SP)
    f5sel <- selectInd(pop = f5, nInd = nSelF6, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    # setting the heritability of the trait, narrow and broad sense
    SP$setVarE(h2 = h2[5], H2 = H2[5])
    f6 <- self(pop = f5sel, nProgeny = 1, simParam = SP)
    f6sel <- selectInd(pop = f6, nInd = 1, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    # setting the heritability of the trait, narrow and broad sense
    SP$setVarE(h2 = h2[6], H2 = H2[6])
    variety <- self(pop = f6sel, nProgeny = 1, simParam = SP)
    
    # updating the TS and marker effects
    TSi <- f3
    popList <- c(popList, TSi)
    if(length(popList) <= 3) {
      TS <- mergePops(popList)
    }
    if(length(popList) > 3) {
      TS <- mergePops(popList[(length(popList)-2):length(popList)]) 
    }
    snps <- RRBLUP2(pop = TS, traits = 1, use = "pheno", snpChip = 1, simParam = SP, maxIter = 20)
    
    newparents <- f4
    pop.base <- newparents
    
    out.Trad_GS <- rbind(out.Trad_GS, data.frame(
      method = "Trad+GS",
      rep = k,
      cycle = i,
      PM = as.numeric(meanG(pop.base)),
      Va = as.numeric(varA(pop.base)),
      Ac = Ac.GS,
      Variety = as.numeric(gv(variety)),
      Years = 5
    ))
  }
}

###########################  saving the outputs  ##############################
cat("saving the output", "\n")
output <- rbind(resultsC0, out.Trad, out.Trad_GS)
write.table(output, "output_simulation.txt")

################################# graphs #######################################
library("tidyr")
library('dplyr')
library("ggplot2")
library('ggpubr')
rm(list = ls()) # clear the environment

data <- read.table("output_simulation.txt")
unique(data$method)
data$method <- as.factor(data$method)
str(data)

PM <-
  data %>%
  group_by(method, cycle) %>% 
  summarise(PM = mean(PM)) %>% 
  ggplot(data = ., mapping = aes(x = cycle, y = PM, col = method)) + 
  geom_line(size = 0.4) +
  geom_point(size = 0.7) +
  labs(x = 'cycles',
       y = 'Population mean') +
  theme(axis.title = element_text(size = 5, face = 'bold'),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 5),
        strip.text = element_text(size = 6, face = 'bold', margin = margin()),
        legend.position = 'bottom',
        legend.text = element_text(size = 5, hjust = 0),
        legend.title = element_text(size = 6, face = 'bold'),
        legend.key.height = unit(3, 'mm'))
PM

##### RS after the breeding timeline - 15 years - 3 breeding cycles
data15 <- data[data$cycle == max(data$cycle),]

p <- ggplot(data15, aes(x = reorder(method, PM), y = PM)) + 
  geom_boxplot(notch = TRUE, outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  geom_jitter(position=position_jitter(0.2)) +
  labs(x = "Breeding method", y = "Response to selection after 15 years of breeding") +
  labs(x = "", fill = "Breeding Method")

p

############################# the end #########################################
attached_packages <- sessionInfo()$otherPkgs
names(attached_packages)