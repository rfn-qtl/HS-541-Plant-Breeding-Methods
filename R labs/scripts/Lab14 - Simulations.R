###################################################
# CS/HS 541 - Plant breeding methods
# Lab 15 - Simulations
# Roberto Fritsche-Neto | roberto.neto@ncsu.edu
# Latest update: June 2, 2026
###################################################

# Load required libraries for simulation, data manipulation, and plotting
library(AlphaSimR)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)

# Record the initial time to track pipeline performance
time.init <- Sys.time()

#################### Parameters ##############
nQTL <- 350                  # Number of Quantitative Trait Loci (QTLs)
add <- 0.22                  # Target mean dominance degree
chip.size.1 <- 910           # SNP chip marker size
segSites <- chip.size.1 + nQTL
replicates <- 5              # Number of simulation replicates (set to 100 for robust research results)
r <- 1:replicates
cycles <- 5                  # Number of breeding cycles per replicate
nChr <- 7                    # Number of chromosomes

# Heritabilities by breeding stage: F2, F3, F4, F5, F6, F7
h2 <- c(0.03, 0.15, 0.40, 0.60, 0.70, 0.80) # Narrow-sense heritability
H2 <- c(0.06, 0.20, 0.45, 0.63, 0.72, 0.81) # Broad-sense heritability

# Crossing block dimensions
n.parents <- 40
n.crosses <- 40
progenie.size <- 100

###########################
# Crop History of Evolution
###########################
set.seed(123)

# Simulate founder historical haplotypes using the Markovian Coalescent Simulator (Macs)
history <- runMacs(nInd = 400, 
                   nChr = nChr,
                   segSites = segSites / nChr, 
                   inbred = TRUE,
                   species = "GENERIC", 
                   split = NULL, 
                   ploidy = 2L)

# Initialize global simulation parameters object
SP = SimParam$new(history)

# Restrict segregating sites to separate QTLs from SNP markers
SP$restrSegSites(nQTL/nChr, chip.size.1/nChr)

# Add an additive + dominance trait with a Gamma distribution for QTL effects
SP$addTraitAD(nQtlPerChr = nQTL / nChr, 
              gamma = TRUE, 
              mean = 0, 
              var = 1, 
              meanDD = add, 
              varDD = 0.5)

# Add a SNP chip layout to the genome
SP$addSnpChip(chip.size.1/nChr)

# Establish environmental variance based on F6/F7 generation baseline heritability
SP$setVarE(h2 = h2[5], H2 = H2[5])

# Generate the initial base population as Doubled Haploids (DH)
firstPop <- newPop(rawPop = history, isDH = TRUE, simParam = SP)

# Check the correlation between Phenotype and True Genetic Value (Accuracy check)
cor(firstPop@pheno, firstPop@gv)

#########################################################
# Simulating the Founders
#########################################################
set.seed(123)

# Phenotypic selection of the top performing individuals to act as initial parents
founders <- selectInd(firstPop, nInd = n.parents, trait = 1, use = "pheno", gender = "B",
                      selectTop = TRUE, returnPop = TRUE, candidates = NULL,
                      simParam = SP)

#############################################################
# Simulating 3 Generations of Traditional Breeding = Burn-In
#############################################################
cat("Running Burn-In Phase...\n")

# Define selection intensities across filial generations
nSelF3 <- n.crosses * progenie.size * 0.1 # Keep top 10%
nSelF4 <- nSelF3 * 0.1
nSelF5 <- nSelF4 * 0.1
nSelF6 <- 1                               # Ultimate variety selection threshold

pop.trad <- founders
perform_list <- list()                    # Using lists instead of rbind for speed

for (i in 1:3) {
  cat("Processing Burn-in Cycle:", i, "\n")
  
  # Year 1: Make random crosses between selected parents
  f1 <- randCross(pop = pop.trad, nCrosses = n.crosses, nProgeny = 1, simParam = SP)
  
  # Year 2: Self-pollinate F1 to create F2 population; perform single plant phenotypic selection
  SP$setVarE(h2 = h2[1], H2 = H2[1])
  f2 <- self(pop = f1, nProgeny = progenie.size, simParam = SP)
  f2sel <- selectInd(pop = f2, nInd = nSelF3, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
  
  # Year 3: Progeny Row Stage 1 (F3 Evaluation)
  SP$setVarE(h2 = h2[2], H2 = H2[2])
  f3 <- self(pop = f2sel, nProgeny = 1, simParam = SP)
  f3sel <- selectInd(pop = f3, nInd = nSelF4, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
  
  # Year 4: Preliminary Yield Test (F4 Evaluation)
  SP$setVarE(h2 = h2[3], H2 = H2[3])
  f4 <- self(pop = f3sel, nProgeny = 1, simParam = SP)
  f4sel <- selectInd(pop = f4, nInd = nSelF5, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
  
  # Year 5: Advanced Test (F5 Evaluation)
  SP$setVarE(h2 = h2[4], H2 = H2[4])
  f5 <- self(pop = f4sel, nProgeny = 1, simParam = SP)
  f5sel <- selectInd(pop = f5, nInd = nSelF6, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
  
  # Years 6-7: Advanced Yield Trials (F6 Evaluation)
  SP$setVarE(h2 = h2[5], H2 = H2[5])
  f6 <- self(pop = f5sel, nProgeny = 1, simParam = SP)
  f6sel <- selectInd(pop = f6, nInd = 1, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
  
  # Year 8: Foundation Seed Increase (F7 Stage)
  SP$setVarE(h2 = h2[6], H2 = H2[6])
  variety <- f6sel
  
  # Recycle parents: F4 generation acts as new parental candidates for the next cycle
  newparents <- f4
  pop.trad <- newparents
  
  # Document genetic progression
  perform_list[[i]] <- data.frame(population = as.numeric(meanG(pop.trad)), 
                                  cultivar = as.numeric(max(variety@gv)))
}

perform <- do.call(rbind, perform_list)
print(perform)

#########################################################
# Creating the initial Training Population (TS0) for GS
#########################################################
cat("Initializing Genomic Selection Training Set (TS0)...\n")

SP$setVarE(h2 = h2[4], H2 = H2[4])
TS0 <- randCross(pop = pop.trad, nCrosses = 24, nProgeny = 48, simParam = SP)
TS00 <- newPop(rawPop = TS0, simParam = SP)

# Generate completely homozygous lines using Doubled Haploid technology
TS000 <- makeDH(pop = TS00, nDH = 1, useFemale = TRUE, simParam = SP)

# Train the genomic prediction model using Ridge Regression BLUP (RR-BLUP)
markers <- RRBLUP2(pop = TS000, traits = 1, use = "pheno", snpChip = 1, simParam = SP, maxIter = 10)

# Evaluate Initial Accuracy of Genomic Selection (GS) on an independent dataset
Testing_set <- randCross(pop = pop.trad, nCrosses = 24, nProgeny = 48, simParam = SP)
Testing_set <- setEBV(Testing_set, markers, simParam = SP)
Ac.GS <- as.numeric(cor(gv(Testing_set), ebv(Testing_set)))
writeLines(paste("Initial GS Accuracy:", Ac.GS), "GS.accuracies.txt")

# Baseline Results DataFrame at Cycle 0
resultsC0 <- data.frame(
  method = c("Trad", "Trad+GS"),
  rep = 1,
  cycle = 0,
  PM = as.numeric(meanG(pop.trad)),
  Va = as.numeric(varA(pop.trad)),
  Ac = c(as.numeric(sqrt(varG(pop.trad)/varP(pop.trad))), Ac.GS),
  Variety = perform[3, 2],
  Years = c(5, 4)
)

print(resultsC0)

##############################################################
# Simulating Breeding Cycles: Traditional Method
##############################################################
cat("Simulating Traditional Breeding Stream...\n")

out.Trad_list <- list()

for (k in 1:length(r)) {
  cat("  -> Replicate:", k, "\n")
  pop.base <- pop.trad
  
  for (i in 1:cycles) {
    f1 <- randCross(pop = pop.base, nCrosses = n.crosses, nProgeny = 1, simParam = SP)
    
    SP$setVarE(h2 = h2[1], H2 = H2[1])
    f2 <- self(pop = f1, nProgeny = progenie.size, simParam = SP)
    f2sel <- selectInd(pop = f2, nInd = nSelF3, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    SP$setVarE(h2 = h2[2], H2 = H2[2])
    f3 <- self(pop = f2sel, nProgeny = 1, simParam = SP)
    f3sel <- selectInd(pop = f3, nInd = nSelF4, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    SP$setVarE(h2 = h2[3], H2 = H2[3])
    f4 <- self(pop = f3sel, nProgeny = 1, simParam = SP)
    f4sel <- selectInd(pop = f4, nInd = nSelF5, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    SP$setVarE(h2 = h2[4], H2 = H2[4])
    f5 <- self(pop = f4sel, nProgeny = 1, simParam = SP)
    f5sel <- selectInd(pop = f5, nInd = nSelF6, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    SP$setVarE(h2 = h2[5], H2 = H2[5])
    f6 <- self(pop = f5sel, nProgeny = 1, simParam = SP)
    f6sel <- selectInd(pop = f6, nInd = 1, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    SP$setVarE(h2 = h2[6], H2 = H2[6])
    variety <- f6sel
    pop.base <- f4 # F4 recycled as new parental population base
    
    idx <- (k - 1) * cycles + i
    out.Trad_list[[idx]] <- data.frame(
      method = "Trad", rep = k, cycle = i,
      PM = as.numeric(meanG(pop.base)),
      Va = as.numeric(varA(pop.base)),
      Ac = as.numeric(sqrt(varG(pop.base)/varP(pop.base))),
      Variety = as.numeric(gv(variety)),
      Years = 5
    )
  }
}
out.Trad <- do.call(rbind, out.Trad_list)

###################################################################
# Simulating Breeding Cycles: Traditional Method + GS
###################################################################
cat("Simulating Traditional + Genomic Selection Stream...\n")

out.Trad_GS_list <- list()

for (k in 1:length(r)) {
  cat("  -> Replicate:", k, "\n")
  
  pop.base <- pop.trad
  TS <- TS000
  popList <- list(TS)
  snps <- markers
  
  for (i in 1:cycles) {
    f1 <- randCross(pop = pop.base, nCrosses = n.crosses, nProgeny = 1, simParam = SP)
    
    SP$setVarE(h2 = h2[1], H2 = H2[1])
    f2 <- self(pop = f1, nProgeny = progenie.size, simParam = SP)
    f2sel <- selectInd(pop = f2, nInd = nSelF3, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    SP$setVarE(h2 = h2[2], H2 = H2[2])
    f3 <- self(pop = f2sel, nProgeny = 1, simParam = SP)
    
    # Genomic Selection implementation at the F3 stage (bypassing phenotypic test wait-times)
    f3.gs <- setEBV(pop = f3, solution = snps, value = "bv", simParam = SP)
    f3sel <- selectInd(pop = f3.gs, nInd = nSelF4, trait = 1, "ebv", gender = "B", selectTop = TRUE, simParam = SP)
    Ac.GS <- as.numeric(cor(gv(f3.gs), ebv(f3.gs)))
    
    SP$setVarE(h2 = h2[3], H2 = H2[3])
    f4 <- self(pop = f3sel, nProgeny = 1, simParam = SP)
    f4sel <- selectInd(pop = f4, nInd = nSelF5, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    SP$setVarE(h2 = h2[4], H2 = H2[4])
    f5 <- self(pop = f4sel, nProgeny = 1, simParam = SP)
    f5sel <- selectInd(pop = f5, nInd = nSelF6, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    SP$setVarE(h2 = h2[5], H2 = H2[5])
    f6 <- self(pop = f5sel, nProgeny = 1, simParam = SP)
    f6sel <- selectInd(pop = f6, nInd = 1, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    SP$setVarE(h2 = h2[6], H2 = H2[6])
    variety <- f6sel
    
    # Update Training Population (TS) with sliding historical window (max 3 cohorts)
    TSi <- f3
    popList <- c(popList, TSi)
    if(length(popList) <= 3) {
      TS <- mergePops(popList)
    } else {
      TS <- mergePops(popList[(length(popList)-2):length(popList)]) 
    }
    snps <- RRBLUP2(pop = TS, traits = 1, use = "pheno", snpChip = 1, simParam = SP, maxIter = 10)
    
    pop.base <- f4
    
    idx <- (k - 1) * cycles + i
    out.Trad_GS_list[[idx]] <- data.frame(
      method = "Trad+GS", rep = k, cycle = i,
      PM = as.numeric(meanG(pop.base)),
      Va = as.numeric(varA(pop.base)),
      Ac = Ac.GS,
      Variety = as.numeric(gv(variety)),
      Years = 4
    )
  }
}
out.Trad_GS <- do.call(rbind, out.Trad_GS_list)

########################### Saving Outputs ##############################
cat("Saving all outputs to text file...\n")
output <- rbind(resultsC0, out.Trad, out.Trad_GS)

# Fixed bug: referenced output$method instead of undefined object method$method
print(table(output$method)) 
write.table(output, "output_simulation.txt", row.names = FALSE)

################################# Graphs #######################################
# Dynamic Clean-up: Remove temporary loops objects but KEEP final output for graphing
rm(list = setdiff(ls(), "output")) 

data <- read.table("output_simulation.txt", header = TRUE)
data$method <- as.factor(data$method)

# Calculate temporal breeding timeline scale
data$Year <- data$cycle * data$Years

# 1. Plotting Genetic Trends (Population Mean over Years)
PM <- data %>%
  group_by(method, Year) %>% 
  summarise(PM = mean(PM, na.rm = TRUE), .groups = "drop") %>% 
  ggplot(aes(x = Year, y = PM, color = method, shape = method, linetype = method, group = method)) + 
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  coord_cartesian(xlim = c(0, 20)) +
  scale_x_continuous(breaks = seq(0, 20, by = 2)) +
  scale_shape_manual(values = c(16, 17)) +   
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(x = "Year", y = "Population Mean", color = "Method", shape = "Method", linetype = "Method") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 8, face = "bold"),
    axis.text = element_text(size = 7),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(PM)

# 2. Response to Selection Boxplot Comparison after a fixed 20-year timeline
data20 <- data %>% filter(Year == 20)

p <- ggplot(data20, aes(x = reorder(method, PM), y = PM, fill = method)) + 
  geom_boxplot(notch = FALSE, outlier.colour = "red", outlier.shape = 16, outlier.size = 2, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 1.5, alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
  scale_fill_manual(values = c("gray70", "gray30")) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  labs(x = "Breeding Methodology", y = "Response to Selection (Year 20)", fill = "Breeding Method") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 9, face = "bold"),
    axis.text = element_text(size = 8),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

print(p)

######## the end ######