## R Version 4.3.2
## 01_imputation.R

# Packages ####
library(stringr) # Version 1.5.1
library(dplyr) # Version 1.1.4
library(rdiversity) # Version 2.2.0
library(ecodist) # Version 2.0.9
library(vegan) # Version 2.6-4
library(doParallel) # Version 1.0.17
library(missForest) # Version 1.5
library(reshape2) # Version 1.4.4
library(sciplot) # Version 1.2-0

# Load and clean data ####
# (Download the BET data from EnviDat: doi.org/10.16904/envidat.348)
BETload <- read.table("./betdata.txt", header=T) 
BETi <- BETload

# Create binary "threatened" variable
unique(BETi$category)
BETi$threatened <- NA
for (n in 1:nrow(BETi)){
  if (isTRUE(BETi$category[n]=="EX"|BETi$category[n]=="RE"|BETi$category[n]=="CR"|BETi$category[n]=="EN"|BETi$category[n]=="VU")){
    BETi$threatened[n] <- "Yes"
  } else if (isTRUE(BETi$category[n]=="NT"|BETi$category[n]=="LC")){
    BETi$threatened[n] <- "No"
  }
}
BETi[c("category","threatened")] <- lapply(BETi[c("category","threatened")], as.factor) 

# Clean data
BETi <- BETi[which(!is.na(BETi$category)==T & BETi$category!="DD"),] # Exclude NA and DD category
BETi[c("EOO_unc","AOO_unc")] <- list(NULL) # Exclude variables EOO_unc and AOO_unc
BETi <- BETi[c(1:9,98,10:97)] # Reorder dataframe

# Recode all traits to correct variable types
# (Excluding taxonomic and Red List columns)
vars <- names(BETi)[11:ncol(BETi)]
# continuous (numeric):
vars_num <- c("size","smind","smeand","smaxd","setmaxl","vpmaxs","lim_low","lim_up","lim_range","EOO_est","AOO_est",
              "nc","MAT","T_diurR","T_iso","T_seas","Tmax_warmM","Tmin_coldM","T_annualR","T_wetQ","T_dryQ","T_warmQ","T_coldQ",
              "MAP","P_wetM","P_dryM","P_seas","P_wetQ","P_dryQ","P_warmQ","P_coldQ","gdd0","gdd5","gdd10","ngd0","ngd5","ngd10")
# ordinal classes:
vars_ord <- c("genl","sfreq","sseas","vpfreq","indL","indT","indK","indF","indR","indN","indS","indHM","biome","eastlim",
              "sub_sum","hab_sum","forest","hemeroby","hem_e")
# nominal categorical:
vars_nom <- c("lstrat","lstrat_e","lform","gform","rK","dmnt","pprot","rhiz","sex","capspos","peri","dwarfm",
              "vp","vp_gem","vp_tub","vp_bul","vp_lea","vp_bra","sub_so","sub_ro","sub_ba","sub_wo","sub_nw","sub_an","epiphyte",
              "aquatic","hab_we","hab_fo","hab_sh","hab_gr","hab_ro","hab_ar")
# Check if all variables are included
sum(length(c(vars_num,vars_ord,vars_nom)))==length(vars) #TRUE
# Apply
BETi[vars_num] <- lapply(BETi[vars_num], as.numeric)
BETi[c(vars_ord,vars_nom)] <- lapply(BETi[c(vars_ord,vars_nom)], as.character)
BETi[c(vars_ord,vars_nom)] <- lapply(BETi[c(vars_ord,vars_nom)], as.factor)

# Overview of missing values
df <- data.frame(
  "Values"=sapply(BETi[13:ncol(BETi)], function(x) length(which(!is.na(x)))),
  "NAs"=sapply(BETi[13:ncol(BETi)], function(x) length(which(is.na(x))))
)

# Exclude variables with >50% missing values
miss50 <- row.names(df[which(df$NAs > df$Values),])
BETi <- BETi[,colnames(BETi)[! colnames(BETi) %in% miss50]]

# Taxonomic eigenvectors ####
tax <- BETi[c(1:6)] # Select taxonomic columns
tax_distance <- c(friendly_name=0, genus=1, family=2, order=3, class=4, phylum=5) # Assign values for each level

# Generate pairwise distances
distm <- rdiversity::tax2dist(tax, tax_distance) # Note: this may take some time
distm2 <- as.matrix(distm@distance)
distm3 <- as.dist(distm2) 
distm3low <- ecodist::lower(distm3)  # Convert symmetric distance matrix to lower triangular matrix

# PCoA: Principal Coordinate Analysis #
set.seed(21)  # Set seed for reproducibility
res <- pco(distm3low)  # Perform Principal Coordinate Analysis
eig <- vegan::eigenvals(res)  # Extract eigenvalues
cum <- cumsum(eig/sum(eig))  # Calculate cumulative proportion explained by each vector
cum[1:20] # Cumulative proportion explained for first 20 eigenvectors 

# Store eigenvectors in dataframe
TVR7 <- as.data.frame(res$vectors[,1:7])
TVR10 <- as.data.frame(res$vectors[,1:10])
TVR15 <- as.data.frame(res$vectors[,1:15])

# Combine trait data with eigenvectors
BETi_TVR7 <- cbind(BETi[c(11:ncol(BETi))], TVR7)
BETi_TVR10 <- cbind(BETi[c(11:ncol(BETi))], TVR10)
BETi_TVR15 <- cbind(BETi[c(11:ncol(BETi))], TVR15)

# Imputation ####
doParallel::registerDoParallel(cores = 8)  # Set up parallel computing with 8 cores

# With 7 eigenvectors
set.seed(91) # Set seed for reproducibility
imp_TVR7 <- missForest(BETi_TVR7, variablewise=F, verbose=T, parallelize = "variables") 

# With 10 eigenvectors
set.seed(99) # Set seed for reproducibility
imp_TVR10 <- missForest(BETi_TVR10, variablewise=F, verbose=T, parallelize = "variables") 

# With 15 eigenvectors
set.seed(99) # Set seed for reproducibility
imp_TVR15 <- missForest(BETi_TVR15, variablewise=F, verbose=T, parallelize = "variables") 

# Without eigenvectors
set.seed(93) # Set seed for reproducibility
imp_TVR0 <- missForest(BETi[c(11:ncol(BETi))], variablewise=F, verbose=T, parallelize="variables")

## Save Out-of-bag errors in dataframe
imp_OOB <- data.frame(imp_TVR0$OOBerror)
imp_OOB <- cbind(
  imp_OOB, imp_TVR7$OOBerror, imp_TVR10$OOBerror, imp_TVR15$OOBerror
)
names(imp_OOB)[1:4] <- c("OOB_TVR0","OOB_TVR7","OOB_TVR10","OOB_TVR15")
imp_OOB <- data.frame(t(imp_OOB))
imp_OOB$TVR <- c(0,7,10,15)
imp_OOB2 <- reshape2::melt(imp_OOB, id="TVR")

# Plot OOB errors
par(mfrow=c(1,2), mar=c(5,5,3,1))
sciplot::lineplot.CI(TVR, value, data=imp_OOB2[which(imp_OOB2$variable=="NRMSE"),], 
            las=1, ylim=c(0.28,0.30), ylab="", xlab="No. of TVRs")
title(ylab = "NRMSE", line=4)
sciplot::lineplot.CI(TVR, value, data=imp_OOB2[which(imp_OOB2$variable=="PFC"),], 
            las=1, ylim=c(0.19,0.21), ylab="", xlab="No. of TVRs")
title(ylab = "PFC", line=4)

# Save imputed data ####
# Select imputation with 7 eigenvectors, save in dataframe
BETi_cmp <- imp_TVR7$ximp
BETi_cmp <- cbind(BETi["threatened"], BETi_cmp[1:ncol(BETi_cmp)]) # Add "threatened" response variable

# For reproducable results: save as .txt file
write.table(BETi_cmp, "BETi_cmp.txt", row.names=F)
#save(BETi_cmp, file="./BETi_cmp.RData")

# End