## R Version 4.3.2
## 02_randomforest.R

## (Run 01_imputation.R first)

# Packages ####
library(Boruta) # Version 8.0.0
library(dplyr) # Version 1.1.4
library(Hmisc) # Version 5.1-1
library(corrplot) # Version 0.92
library(DescTools) # Version 0.99.52
library(randomForest) # Version 4.7-1.1
library(rfUtilities) # Version 2.1-5

# Load and clean data ####
# Load imputed dataset from 01_imputation.R
# For complete reproducability, use .txt file
BETi_cmp <- read.table("./BETi_cmp.txt", header=T)
#load(file="BETi_cmp.RData")

# Save in new dataframe
BETi_cmp2 <- BETi_cmp

# Recode to correct variable types
vars_num2 <- vars_num[! vars_num %in% miss50]
BETi_cmp2[vars_num2] <- lapply(BETi_cmp2[vars_num2], as.numeric)
vars_ord2 <- vars_ord[! vars_ord %in% miss50]
vars_nom2 <- vars_nom[! vars_nom %in% miss50]
BETi_cmp2[c(vars_ord2,vars_nom2)] <- lapply(BETi_cmp2[c(vars_ord2,vars_nom2)], as.character)
BETi_cmp2[c(vars_ord2,vars_nom2)] <- lapply(BETi_cmp2[c(vars_ord2,vars_nom2)], as.factor)
BETi_cmp2$threatened <- as.factor(BETi_cmp2$threatened)

# Exclude variables related to range size
BETi_cmp2[c("EOO_est","AOO_est","nc")] <- list(NULL) 

# Organize into different variable groups
names(BETi_cmp2)
vars_bio <- c("size","lstrat","lstrat_e","lform","gform","genl","rK","dmnt",
              "sex","smind","smeand","smaxd","sfreq","sseas","capspos","peri",
              "vp","vp_gem","vp_tub","vp_bul","vp_lea","vp_bra")
vars_eco <- c("indT","indK","biome","eastlim","lim_low","lim_up","lim_range",
              "indL","indF","indR","indN","indS","indHM",
              "sub_so","sub_ro","sub_ba","sub_wo","sub_nw","sub_an","sub_sum","epiphyte","aquatic",
              "hab_we","hab_fo","hab_sh","hab_gr","hab_ro","hab_ar","hab_sum","hem_e")
vars_cli <- c("MAT","T_diurR","T_iso","T_seas","Tmax_warmM","Tmin_coldM","T_annualR","T_wetQ","T_dryQ","T_warmQ","T_coldQ",
              "MAP","P_wetM","P_dryM","P_seas","P_wetQ","P_dryQ","P_warmQ","P_coldQ",
              "gdd0","gdd5","gdd10","ngd0","ngd5","ngd10")
vars_TVR <- c("V1","V2","V3","V4","V5","V6","V7")
# Check if all variables are included
length(c("threatened", vars_bio, vars_eco, vars_cli, vars_TVR)) == ncol(BETi_cmp2) #TRUE

# Store in separate dataframes
tt_all <- BETi_cmp2[c("threatened", vars_bio, vars_eco, vars_cli, vars_TVR)]
tt_bio <- BETi_cmp2[c("threatened", vars_bio, vars_TVR)]
tt_eco <- BETi_cmp2[c("threatened", vars_eco, vars_TVR)]
tt_cli <- BETi_cmp2[c("threatened", vars_cli, vars_TVR)]

# Variable selection and Random Forest ####
## Combined (all trait groups) ####
# Feature selection with Boruta
set.seed(12)
ba <- Boruta::Boruta(threatened ~., data=tt_all)
vars_bor <- Boruta::getSelectedAttributes(ba, withTentative=TRUE) # Select variables that are confirmed or tentative
tt_all2 <- tt_all[c("threatened", vars_bor)] # Keep selected variables

# Divide into continuous and categorical variables
num <- names(tt_all2[2:ncol(tt_all2)] %>% dplyr::select(where(is.numeric)))
fac <- names(tt_all2[2:ncol(tt_all2)] %>% dplyr::select(where(is.factor)))
tt_all2n <- tt_all2[num]
tt_all2f <- tt_all2[fac]

# Check for correlations > 0.7
# Continuous variables: Pearson
set.seed(123)
rc <- rcorr(as.matrix(tt_all2n), type="pearson")
par(mfrow = c(1,1), mar = rep(0.1, 4))
corrplot.mixed(rc$r, tl.cex=0.5, number.cex=0.5, tl.col="black")

# Categorical variables: Cramer's V
set.seed(123)
cor.v <- PairApply(tt_all2f, CramerV)
par(mfrow = c(1,1), mar = rep(0.1, 4))
corrplot.mixed(cor.v, tl.cex=0.5, number.cex=0.5, tl.col="black")

# Select correlated variables
vars_cor <- c("lstrat_e","sub_ba","sub_wo","smind","smaxd","vp_gem","genl","rK","peri",
              "lim_up","gdd0","gdd5","gdd10","ngd0","ngd5","ngd10",
              "Tmin_coldM","Tmax_warmM","T_dryQ","T_warmQ","T_coldQ",
              "T_annualR","P_wetM","P_dryM","P_wetQ","P_dryQ","P_coldQ")
vars_sl <- colnames(tt_all2)[! colnames(tt_all2) %in% vars_cor]
# Exclude correlated variables
tt_all3 <- tt_all2[vars_sl]

# Random Forest modeling
# Change order of factor levels for "threatened"
tt_all3$threatened <- factor(tt_all3$threatened, levels=c("Yes","No")) 

# Apply weights to balance out classes
data_wt <- ifelse(tt_all3$threatened=="Yes", 1/(389/nrow(tt_all3)), 1/((nrow(tt_all3)-389)/nrow(tt_all3))) #weights

# Construct Random Forest 
set.seed(51)
rf.all <- randomForest::randomForest(
  tt_all3[-1], tt_all3[[1]], weights=data_wt, 
  importance=TRUE, keep.forest=TRUE,
  ntree=500, mtry=floor(sqrt(ncol(tt_all3[-1]))) 
  ) 
# Show model
rf.all

# Get True Positive, False Positive, True Negative and False Negative
table(predict(rf.all))
TP <- rf.all$confusion[1] 
FP <- rf.all$confusion[2] 
FN <- rf.all$confusion[3] 
TN <- rf.all$confusion[4] 

# Calculate performance metrics
rfUtilities::accuracy(cbind(c(TP,FP),c(FN,TN))) 
all_pre <- TP/(TP+FP) # Precision

## Biological traits ####
# Feature selection with Boruta
set.seed(12)
bb <- Boruta::Boruta(threatened ~., data=tt_bio) #threatened
vars_bor <- Boruta::getSelectedAttributes(bb, withTentative=TRUE) # Select variables that are confirmed or tentative
tt_bio2 <- tt_bio[c("threatened", vars_bor)] # Keep selected variables

# Divide into continuous and categorical variables
num <- names(tt_bio2[2:ncol(tt_bio2)] %>% select(where(is.numeric)))
fac <- names(tt_bio2[2:ncol(tt_bio2)] %>% select(where(is.factor)))
tt_bio2n <- tt_bio2[num]
tt_bio2f <- tt_bio2[fac]

# Check for correlations > 0.7
# Continuous variables: Pearson
set.seed(123)
rc <- rcorr(as.matrix(tt_bio2n), type="pearson")
par(mfrow = c(1,1), mar = rep(0.1, 4))
corrplot.mixed(rc$r, tl.cex=0.5, number.cex=0.5, tl.col="black") #remove: smind, smaxd

# Categorical variables: Cramer's V
set.seed(123)
cor.v <- PairApply(tt_bio2f, CramerV)
par(mfrow = c(1,1), mar = rep(0.1, 4))
corrplot.mixed(cor.v, tl.cex=0.5, number.cex=0.5, tl.col="black") #remove: lstrat_e, genl, rK, peri

# Select correlated variables
vars_cor <- c("smind","smaxd","lstrat_e","vp_gem","genl","rK","peri")
vars_sl <- colnames(tt_bio2)[! colnames(tt_bio2) %in% vars_cor]
# Exclude correlated variables
tt_bio3 <- tt_bio2[vars_sl]

# Random Forest modeling
# Change order of factor levels for "threatened"
tt_bio3$threatened <- factor(tt_bio3$threatened, levels=c("Yes","No")) #change order of threatened levels

# Apply weights to balance out classes
data_wt <- ifelse(tt_bio3$threatened=="Yes", 1/(389/nrow(tt_bio3)), 1/((nrow(tt_bio3)-389)/nrow(tt_bio3))) #giving threatened species larger weights (balanced threatened/non-threatened)

# Construct Random Forest
set.seed(51)
rf.bio <- randomForest::randomForest(
  tt_bio3[-1], tt_bio3[[1]], weights=data_wt, 
  importance=TRUE, keep.forest=TRUE,
  ntree=500, mtry=floor(sqrt(ncol(tt_bio3[-1]))) 
  )
# Show model
rf.bio 

# Get True Positive, False Positive, True Negative and False Negative
table(predict(rf.bio))
TP <- rf.bio$confusion[1] 
FP <- rf.bio$confusion[2] 
FN <- rf.bio$confusion[3] 
TN <- rf.bio$confusion[4] 

# Calculate performance metrics
rfUtilities::accuracy(cbind(c(TP,FP),c(FN,TN))) 
bio_acc <- (TP+TN)/(TP+TN+FP+FN) # Accuracy
bio_spe <- TN/(TN+FP) # Specificity
bio_sen <- TP/(TP+FN) # Sensitivity
bio_pre <- TP/(TP+FP) # Precision
bio_fpr <- FP/(FP+TN) # False Positive Rate
bio_auc <- (bio_sen - bio_fpr +1)/2 #AUC

## Ecological traits ####
# Feature selection with Boruta
set.seed(12)
be <- Boruta::Boruta(threatened ~., data=tt_eco) #threatened
vars_bor <- Boruta::getSelectedAttributes(be, withTentative=TRUE) # Select variables that are confirmed or tentative
tt_eco2 <- tt_eco[c("threatened", vars_bor)] # Keep selected variables

# Divide into continuous and categorical variables
num <- names(tt_eco2[2:ncol(tt_eco2)] %>% select(where(is.numeric)))
fac <- names(tt_eco2[2:ncol(tt_eco2)] %>% select(where(is.factor)))
tt_eco2n <- tt_eco2[num]
tt_eco2f <- tt_eco2[fac]

# Check for correlations > 0.7
# Continuous variables: Pearson
set.seed(123)
rc <- rcorr(as.matrix(tt_eco2n), type="pearson")
par(mfrow = c(1,1), mar = rep(0.1, 4))
corrplot.mixed(rc$r, tl.cex=0.5, number.cex=0.5, tl.col="black")

# Categorical variables: Cramer's V
set.seed(123)
cor.v <- PairApply(tt_eco2f, CramerV)
par(mfrow = c(1,1), mar = rep(0.1, 4))
corrplot.mixed(cor.v, tl.cex=0.5, number.cex=0.5, tl.col="black")

# Select correlated variables
vars_cor <- c("lim_up","sub_ba","sub_wo")
vars_sl <- colnames(tt_eco2)[! colnames(tt_eco2) %in% vars_cor]
# Exclude correlated variables
tt_eco3 <- tt_eco2[vars_sl]

# Random Forest
# Change order of factor levels for "threatened"
tt_eco3$threatened <- factor(tt_eco3$threatened, levels=c("Yes","No")) #change order of threatened levels

# Apply weights to balance out classes
data_wt <- ifelse(tt_eco3$threatened=="Yes", 1/(389/nrow(tt_eco3)), 1/((nrow(tt_eco3)-389)/nrow(tt_eco3))) #give threatened species larger weights

# Construct Random Forest
set.seed(52)
rf.eco <- randomForest::randomForest(
  tt_eco3[-1], tt_eco3[[1]], weights=data_wt, 
  importance=TRUE, keep.forest=TRUE,
  ntree=500, mtry=floor(sqrt(ncol(tt_eco3[-1])))
  )
# Show model
rf.eco 

# Get True Positive, False Positive, True Negative and False Negative
table(predict(rf.eco))
TP <- rf.eco$confusion[1] 
FP <- rf.eco$confusion[2] 
FN <- rf.eco$confusion[3] 
TN <- rf.eco$confusion[4] 

# Calculate performance metrics
rfUtilities::accuracy(cbind(c(TP,FP),c(FN,TN))) 
eco_pre <- TP/(TP+FP) # Precision

## Bioclimatic variables ####
# Feature selection with Boruta
set.seed(12)
bc <- Boruta(threatened ~., data=tt_cli) #threatened
vars_bor <- Boruta::getSelectedAttributes(bc, withTentative=TRUE) # Select variables that are confirmed or tentative
tt_cli2 <- tt_cli[c("threatened", vars_bor)] # Keep selected variables

# Check for correlations > 0.7
# Continuous variables: Pearson
set.seed(123)
rc <- rcorr(as.matrix(tt_cli2[2:ncol(tt_cli2)]), type="pearson")
par(mfrow = c(1,1), mar = rep(0.1, 4))
corrplot.mixed(rc$r, tl.cex=0.5, number.cex=0.5, tl.col="black") #

# Select correlated variables
vars_cor <- c("gdd0","gdd5","gdd10","ngd0","ngd5","ngd10",
              "Tmin_coldM","Tmax_warmM","T_dryQ","T_warmQ","T_coldQ",
              "T_annualR", "P_wetM","P_dryM","P_wetQ","P_dryQ","P_coldQ")
vars_sl <- colnames(tt_cli2)[! colnames(tt_cli2) %in% vars_cor]
# Exclude correlated variables
tt_cli3 <- tt_cli2[vars_sl]

# Random Forest 
# Change order of factor levels for "threatened"
tt_cli3$threatened <- factor(tt_cli3$threatened, levels=c("Yes","No")) 

# Apply weights to balance out classes
data_wt <- ifelse(tt_cli3$threatened=="Yes", 1/(389/nrow(tt_cli3)), 1/((nrow(tt_cli3)-389)/nrow(tt_cli3))) 

# Construct Random Forest
set.seed(51)
rf.cli <- randomForest::randomForest(
  tt_cli3[-1], tt_cli3[[1]], weights=data_wt, 
  importance=TRUE, keep.forest=TRUE,
  ntree=500, mtry=floor(sqrt(ncol(tt_cli3[-1]))) 
  )
# Show model
rf.cli

# Get True Positive, False Positive, True Negative and False Negative
table(predict(rf.cli))
TP <- rf.cli$confusion[1] 
FP <- rf.cli$confusion[2] 
FN <- rf.cli$confusion[3] 
TN <- rf.cli$confusion[4] 

# Calculate performance metrics
rfUtilities::accuracy(cbind(c(TP,FP),c(FN,TN))) 
cli_pre <- TP/(TP+FP) # Precision

# End