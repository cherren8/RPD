#### Supplementary Figures for rpd Methods Paper. Begins with script used in main text:

library(vegan)
library(moments)
library(vioplot)
library(pROC)

####
#### Breast health

# This version is for picking OTUs with 100nt and a read cutoff of 10k
b <- read.csv("~/Desktop/BreastHealthFiles/BreastHealthPickOTUcutoff10k_08Apr26.csv", header = T)

# Read in mapping file
m <- read.delim("~/Desktop/BreastHealthFiles/BreastHealthMetadata.txt", header = T)
colnames(m)

# Remove blanks
b <- b[!grepl("BLANK", b$X), ]
dim(b)

# Remove OTU names to leave jut numeric matrix
b.num <- b[, -1]

# Create a function to pull out one item from a vector
select.vec <- function(n, vec){
  back <- vec[n]
  return(back)
}

# Match names of m and b
samp.id1 <- sapply(strsplit(m$sample_name, "[.]"), select.vec, n=1)
samp.id2 <- sapply(strsplit(m$sample_name, "[.]"), select.vec, n=2)
m$names <- paste(samp.id1, samp.id2, sep = ".")

m.ord <- match(b$X, m$names)
length(m.ord)
m.sub <- m[m.ord, ]

# Check that number of rows are equal
dim(m.sub)[1] == dim(b.num)[1]

## Create reference, case, and control datasets

# Count how many cases and controls there are in matching mapping data
table(m.sub$case_control__rtl_case_or_contro)

# Pull out all controls and randomize them
healthy.seed <- 1111 #gave good results for gbh
set.seed(healthy.seed)
healthy.gbh <- sample(which(m.sub$case_control__rtl_case_or_contro == "control"))

# Create reference set from first n controls
ref.size <- 100
ref.threshold <- .7
ref.index <- 1:ref.size
ref <- healthy.gbh[ref.index]
# Create control dataset from remainder of healthy dataset
control <- healthy.gbh[-ref.index]
# Create case dataset from samples labeled "case"
case <- which(m.sub$case_control__rtl_case_or_contro == "case")

# Calculate rpd for controls
control.rpd.gbh <- calc.rpd(ref = b.num[ref, ], test = b.num[control, ], ref.pers = ref.threshold)

# Calculate rpd for cases
case.rpd.gbh <- calc.rpd(ref = b.num[ref, ], test = b.num[case, ], ref.pers = ref.threshold)

case.control.vec.gbh <- c(rep("case", times = length(case)), rep("control", times = length(control)))

table(case.control.vec.gbh)

rpd.vec.gbh <- c(case.rpd.gbh, control.rpd.gbh)

# Calculate summary stats on rpd values
mean(case.rpd.gbh)
sd(case.rpd.gbh)

mean(control.rpd.gbh)
sd(control.rpd.gbh)

# Convert case/control vector to 1/0
case.control.binary.gbh <- as.numeric(case.control.vec.gbh == "case")

# Implement logistic regression
summary(glm(as.factor(case.control.binary.gbh) ~ rpd.vec.gbh, family = "binomial"))

## Calculate relative risk of cases in the top XX% of rpd values
risk.fract <- .5

# Find number of observations this fraction corresponds to
num.order <- round(risk.fract * length(case.control.vec.gbh))

# Find the indices for the top fraction of rpd values
top.fract <- tail(order(rpd.vec.gbh, na.last = NA), n = num.order)

# Calculate risk ratio
mean(case.control.vec.gbh[top.fract] == "case") / mean(case.control.vec.gbh[-top.fract] == "case")

# Check how many OTUs retained at threshold
sum(colMeans(ceiling(b.num[ref, ])) > ref.threshold)
b.num.ref.sub <- b.num[ref, colMeans(ceiling(b.num[ref, ])) > ref.threshold]
mean(rowSums(b.num.ref.sub))

## Repeat analysis using BCD to reference centroid

# Create centroid by averaging reference samples
# bc.ref.mean.gbh <- apply(b.num[ref,], 2, mean)
# # Create matrix of case and control samples for use in Bray-Curtis analysis
# b.num.bc <- rbind(b.num[case, ], b.num[control, ])
#
# # Calculate BCD to centroid for each sample
# bc.vec.ref <- as.matrix(vegdist(rbind(bc.ref.mean.gbh, b.num.bc)))[-1, 1]
#
# # Implement same logistic regression
# summary(glm(as.factor(case.control.binary.gbh) ~ bc.vec.ref, family = "binomial"))
#
# ## Calculate diversity metrics
# ref.div.gbh <- vegan::diversity(b.num[ref, ])
# control.div.gbh <- vegan::diversity(b.num[control, ])
# case.div.gbh <- vegan::diversity(b.num[case, ])
# case.control.div.gbh <- c(case.div.gbh, control.div.gbh)

####
#### American Gut Data

d <- read.csv("~/Desktop/HarvardDSI/Projects/ControlMatching/Data/AGPpreprint/AGPpreprintRename5_20Apr18.csv", header = T, row.names = 1)
dim(d)
meta <- read.csv("~/Desktop/HarvardDSI/Projects/ControlMatching/Data/AGPpreprint/AGPmapping_02May18.csv", header = T, row.names = 1)


#Subset and put metadata into same order as d
meta.ord <- match(rownames(d), rownames(meta))
length(meta.ord)
meta.sub <- meta[meta.ord, ]

#Check that samples are in the same order
mean(rownames(meta.sub) == rownames(d))

#Take out categories from metadata that do not have sufficient variation to be useful as a predictor
#Make a function to indicate whether more than XX% of a vector is the same factor
prop.same <- function(vec) {
  prop.largest <- (max(table(vec)) / length(vec) )
  return(prop.largest)
}

summary.max <- function(vec){
  (max(table(vec)))
}

# Remove the vioscreen diet columns
vio.col <- unique(c(grep("vioscreen", colnames(m.sub)) ) )

# Identify which columns are too homogeneous to use as a predictor (70% or more the same value)
prop.same.col <- which(apply(meta.sub, 2, prop.same) > .7)
out.col <- unique(c(vio.col, prop.same.col))

#Take out people with no year or weight and who are from outside USA, and non-stool samples
out.row <- unique(c(grep("provided", meta.sub$age_years), grep("provided", meta.sub$weight_kg) , grep("provided", meta.sub$bmi) , which(meta.sub$country != "USA") , which(meta.sub$sample_type != "Stool") ) )

# Remove unnecessary columns and matching rows from metadata
meta.sub1 <- meta.sub[-out.row, -out.col]
# Remove rows from stool data
d.sub1 <- d[-out.row, ]

# Check dimension of d.sub1
dim(d.sub1)

# Check again that row names are the same
mean(rownames(d.sub1) == rownames(meta.sub1))

par(mfrow = c(1, 1))

# Create reference dataset of healthy individuals

set.seed(healthy.seed)
healthy <- sample(which(meta.sub1$cdiff == "I do not have this condition"
                        & meta.sub1$ibs == "I do not have this condition"
                        & meta.sub1$antibiotic_history == "I have not taken antibiotics in the past year."
                        & meta.sub1$cancer == "I do not have this condition"
                        & !meta.sub1$age_cat %in% c("child", "teen")) )

length(healthy)
dim(meta.sub1)

# Take first 500 healthy to create a reference
healthy.ref <- healthy[1:500]

# Pull out another set of healthy, plus antibiotic individuals
ab.test.indices <- c(healthy[501:length(healthy)],
                     which(meta.sub1$antibiotic_history == "Week"),
                     which(meta.sub1$antibiotic_history == "Month"))

table(meta.sub1$antibiotic_history)

# Create reference dataset
ab.ref <- as.matrix(d.sub1[healthy.ref, ])
# Create test dataset (contains both controls and cases)
ab.test <- as.matrix(d.sub1[ab.test.indices, ])
# Create labels for test dataset
meta.ab.test <- meta.sub1$antibiotic_history[ab.test.indices]

ab.ref.pers <- .7

# Calculate rpd for antibiotic test dataset
ab.test.rpd <- calc.rpd(ab.ref, ab.test, ref.pers = ab.ref.pers)

# Check how many taxa this cutoff keeps
sum(colMeans(ceiling(ab.ref)) > ab.ref.pers)

ab.ref.sub <- ab.ref[, colMeans(ceiling(ab.ref)) > ab.ref.pers]
mean(rowSums(ab.ref.sub))

# Calculate summary stats on rpds
mean(ab.test.rpd[grepl("not", meta.ab.test)])
sd(ab.test.rpd[grepl("not", meta.ab.test)])

mean(ab.test.rpd[!grepl("not", meta.ab.test)], na.rm = T)
sd(ab.test.rpd[!grepl("not", meta.ab.test)], na.rm = T)

# Make case/control vector into 1/0 vector
case.control.binary.abx <- as.numeric(meta.ab.test == "Week" | meta.ab.test == "Month" )

# Implement logistic regression
summary(glm(as.factor(case.control.binary.abx) ~ ab.test.rpd, family = "binomial"))

## Do Risk Ratio analysis
# Find number of observations this fraction corresponds to
num.order.abx <- round(risk.fract * length(case.control.binary.abx))

# Find the indices for the top fraction of rpd values
top.fract.abx <-  tail(order(ab.test.rpd, na.last = NA), n = num.order.abx)

# Calculate risk ratio
mean(case.control.binary.abx[top.fract.abx] == 1) / mean(case.control.binary.abx[-top.fract.abx] == 1)

## Repeat analysis wit BCD to reference centroid
# Calculate centroid for ref group
# ab.ref.mean <- apply(ab.ref, 2, mean)
# # Calculate BCD to centroid
# bc.vec.ab <- as.matrix(vegdist(rbind(ab.ref.mean, ab.test)))[-1, 1]
#
# # Repeat logistic regression
# summary(glm(as.factor(case.control.binary.abx) ~ bc.vec.ab, family = "binomial"))

#### C diff

# Check for overlap between C. diff and antibiotics groups
table(meta.sub1$cdiff)
sum(meta.sub1$cdiff == "Diagnosed by a medical professional (doctor, physician assistant)")
# Check how many of these people took antibiotics recently
sum(meta.sub1$cdiff == "Diagnosed by a medical professional (doctor, physician assistant)" &
      meta.sub1$antibiotic_history %in% c("Month", "Week") ) # groups largely do not overlap

# Pull out same set of healthy controls, plus cdiff individuals
cd.test.indices <- c(healthy[501:length(healthy)],
                     which(meta.sub1$cdiff == "Diagnosed by a medical professional (doctor, physician assistant)"))

# Create reference and test datasets, and labels for test
cd.ref <- as.matrix(d.sub1[healthy.ref, ])
cd.test <- as.matrix(d.sub1[cd.test.indices, ])
meta.cd.test <- meta.sub1$cdiff[cd.test.indices]

cd.ref.pers <- .7

# Calculate rpd
cd.test.rpd <- calc.rpd(cd.ref, cd.test, ref.pers = cd.ref.pers)

# Check how many taxa are retained
sum(colMeans(ceiling(cd.ref)) > cd.ref.pers)
cd.ref.sub <- cd.ref[, colMeans(ceiling(cd.ref)) > cd.ref.pers]
mean(rowSums(cd.ref.sub))

# Calculate summary stats on cd samples

mean(cd.test.rpd[!grepl("Diagnosed", meta.cd.test)], na.rm = T)
sd(cd.test.rpd[!grepl("Diagnosed", meta.cd.test)], na.rm = T)

mean(cd.test.rpd[grepl("Diagnosed", meta.cd.test)], na.rm = T)
sd(cd.test.rpd[grepl("Diagnosed", meta.cd.test)], na.rm = T)

# Create 0/1 vector
case.control.binary.cd <- as.numeric(grepl("Diagnosed", meta.cd.test ))

# Implement logistic regression
summary(glm(as.factor(case.control.binary.cd) ~ cd.test.rpd, family = "binomial"))

# Do Relative Risk analysis
# Find number of observations this fraction corresponds to
num.order.cd <- round(risk.fract * length(case.control.binary.cd))

# Find the indices for the top fraction of rpd values
top.fract.cd <-  tail(order(cd.test.rpd, na.last = NA), n = num.order.cd)

# Calculate relative risk
mean(case.control.binary.cd[top.fract.cd] == 1) / mean(case.control.binary.cd[-top.fract.cd] == 1)

## Repeat analysis using BCD to centroid
# Calculate centroid
# cd.ref.mean <- apply(cd.ref, 2, mean)
#
# bc.vec.cd <- as.matrix(vegdist(rbind(cd.ref.mean, cd.test)))[-1, 1]
#
# # Implement logistic regression
# summary(glm(as.factor(case.control.binary.cd) ~ bc.vec.cd, family = "binomial"))

####################################################################################################
####################################################################################################

#### AUC depending on ref.pers

## Set seed to make reproducible
set.seed(4444)

#### Set reference size
ref.size <- 50
#### Choose different ref pers thresholds
ref.threshold.vec <- seq(.6, .95, .05)
#### Set number of iterations for each ref.pers value
n.ref.pers.iter <- 50

## GBH

# Begin as before:

# Pull out all controls
healthy.gbh <- which(m.sub$case_control__rtl_case_or_contro == "control")

for(entry in ref.threshold.vec){

  ref.thres <- entry

  #### For loop to randomize reference data
  auc.ref.pers.gbh <- vector()

  for(iter in 1:n.ref.pers.iter){

    # Randomize indices
    ref.index <- sample(1:length(healthy.gbh), ref.size)
    ref <- healthy.gbh[ref.index]
    # Create control dataset from remainder of healthy dataset
    control <- healthy.gbh[-ref.index]
    # Create case dataset from samples labeled "case"
    case <- which(m.sub$case_control__rtl_case_or_contro == "case")

    # Calculate rpd
    rpd.vec.gbh <- calc.rpd(ref = as.matrix(b.num[ref, ]),
                            test = as.matrix(b.num[c(case, control), ]),
                            ref.pers = ref.thres)

    case.control.vec.gbh <- c(rep("case", times = length(case)), rep("control", times = length(control)))

    # Convert case/control vector to 1/0
    case.control.binary.gbh <- as.numeric(case.control.vec.gbh == "case")

    # Implement logistic regression
    gbh.data.ref.pers <- as.data.frame(cbind(case.control.binary.gbh, rpd.vec.gbh))
    model.br.rpd <- glm(case.control.binary.gbh ~ rpd.vec.gbh , data = gbh.data.ref.pers, family = "binomial")

    # Generate predicted probabilities
    preds.br.rpd <- predict(model.br.rpd, newdata = gbh.data.ref.pers, type = "response")

    # Create and plot the ROC curve
    roc_obj.br.rpd <- roc(gbh.data.ref.pers$case.control.binary.gbh, preds.br.rpd)

    auc.ref.pers.gbh[iter] <- auc(roc_obj.br.rpd)

  }

  results.pers.name.gbh <- paste0("auc.gbh.pers.loop", ref.thres)

  assign(results.pers.name.gbh, auc.ref.pers.gbh)

  hist(auc.ref.pers.gbh)

}



# May run faster for AGP since there are fewer taxa at persistence thresholds

####################################################################################################
####################################################################################################

#### AGP Antibiotics

# Create reference dataset of healthy individuals

healthy <- which(meta.sub1$cdiff == "I do not have this condition"
                 & meta.sub1$ibs == "I do not have this condition"
                 & meta.sub1$antibiotic_history == "I have not taken antibiotics in the past year."
                 & meta.sub1$cancer == "I do not have this condition"
                 & !meta.sub1$age_cat %in% c("child", "teen")
)

length(healthy)
dim(meta.sub1)


#### Set reference size
ref.size <- 200 # Larger takes longer to run. Main text analysis used 500
#### Choose different ref pers thresholds
ref.threshold.vec.agp <- seq(.5, .9, .05)
#### Set number of iterations for each ref.pers value
n.ref.pers.iter.agp <- 100

set.seed(4444)

for(entry in ref.threshold.vec.agp){

  ref.thres.agp <- entry


  #### For loop to randomize reference data
  auc.ref.pers.agp.ab <- vector()
  auc.ref.pers.agp.cd <- vector()

  for(iter in 1:n.ref.pers.iter.agp){

    # Randomize indices
    ref.index <- sample(1:length(healthy), ref.size)
    ref <- healthy[ref.index]
    # Create control dataset from remainder of healthy dataset
    control <- healthy[-ref.index]
    # Create case dataset from samples labeled "case"
    case <- which(m.sub$case_control__rtl_case_or_contro == "case")


    # Take sample of healthy to create a reference
    healthy.ref.indices <- sample(1:length(healthy), ref.size)

    healthy.ref <- healthy[healthy.ref.indices]

    healthy.control <- sample(healthy[-healthy.ref.indices], 199) # 199 is # left when ref size is 500, so use this for all

    # Pull out another set of healthy, plus antibiotic individuals
    ab.test.indices <- c(healthy.control,
                         which(meta.sub1$antibiotic_history == "Week"),
                         which(meta.sub1$antibiotic_history == "Month"))

    table(meta.sub1$antibiotic_history)

    # Create reference dataset
    ab.ref <- as.matrix(d.sub1[healthy.ref, ])
    # Create test dataset (contains both controls and cases)
    ab.test <- as.matrix(d.sub1[ab.test.indices, ])
    # Create labels for test dataset
    meta.ab.test <- meta.sub1$antibiotic_history[ab.test.indices]

    # Calculate rpd for antibiotic test dataset
    ab.test.rpd <- calc.rpd(ab.ref, ab.test, ref.pers = ref.thres.agp)

    # Make case/control vector into 1/0 vector
    case.control.binary.abx <- as.numeric(meta.ab.test == "Week" | meta.ab.test == "Month" )

    # Implement logistic regression
    agp.ab.data.ref.pers <- as.data.frame(cbind(case.control.binary.abx, ab.test.rpd))
    model.agp.ab.rpd <- glm(case.control.binary.abx ~ ab.test.rpd , data = agp.ab.data.ref.pers, family = "binomial")

    # Generate predicted probabilities
    preds.agp.ab.rpd <- predict(model.agp.ab.rpd, newdata = agp.ab.data.ref.pers, type = "response")

    # Create and plot the ROC curve
    roc_obj.agp.ab.rpd <- roc(agp.ab.data.ref.pers$case.control.binary.abx, preds.agp.ab.rpd)

    auc.ref.pers.agp.ab[iter] <- auc(roc_obj.agp.ab.rpd)


    ## Repeat for C. diff

    # Pull out another set of healthy, plus antibiotic individuals
    cd.test.indices <- c(healthy.control,
                         which(meta.sub1$cdiff == "Diagnosed by a medical professional (doctor, physician assistant)"))

    # Create reference dataset
    cd.ref <- as.matrix(d.sub1[healthy.ref, ])
    # Create test dataset (contains both controls and cases)
    cd.test <- as.matrix(d.sub1[cd.test.indices, ])
    # Create labels for test dataset
    meta.cd.test <- meta.sub1$cdiff[cd.test.indices]

    # Calculate rpd for antibiotic test dataset
    cd.test.rpd <- calc.rpd(cd.ref, cd.test, ref.pers = ref.thres.agp)

    # Create 0/1 vector
    case.control.binary.cd <- as.numeric(grepl("Diagnosed", meta.cd.test ))

    # Implement logistic regression
    agp.cd.data.ref.pers <- as.data.frame(cbind(case.control.binary.cd, cd.test.rpd))
    model.agp.cd.rpd <- glm(case.control.binary.cd ~ cd.test.rpd , data = agp.cd.data.ref.pers, family = "binomial")

    # Generate predicted probabilities
    preds.agp.cd.rpd <- predict(model.agp.cd.rpd, newdata = agp.cd.data.ref.pers, type = "response")

    # Create and plot the ROC curve
    roc_obj.agp.cd.rpd <- roc(agp.cd.data.ref.pers$case.control.binary.cd, preds.agp.cd.rpd)

    auc.ref.pers.agp.cd[iter] <- auc(roc_obj.agp.cd.rpd)


  }

  par(mfrow = c(1, 2))

  results.pers.name.agp.ab <- paste0("auc.agp.ab.pers.loop", ref.thres.agp)

  assign(results.pers.name.agp.ab, auc.ref.pers.agp.ab)

  hist(auc.ref.pers.agp.ab)

  results.pers.name.agp.cd <- paste0("auc.agp.cd.pers.loop", ref.thres.agp)

  assign(results.pers.name.agp.cd, auc.ref.pers.agp.cd)

  hist(auc.ref.pers.agp.cd)
}

## Note that these AUC values might be lower than the text values because the reference size is 100 in order to reduce overlap between randomly sampled reference sets

par(mfrow = c(1, 1))

mean(auc.gbh.pers.loop0.6)
mean(auc.gbh.pers.loop0.65)
mean(auc.gbh.pers.loop0.7)
mean(auc.gbh.pers.loop0.75)
mean(auc.gbh.pers.loop0.8)
mean(auc.gbh.pers.loop0.85)
mean(auc.gbh.pers.loop0.9)
mean(auc.gbh.pers.loop0.95)

sd(auc.gbh.pers.loop0.6)
sd(auc.gbh.pers.loop0.65)
sd(auc.gbh.pers.loop0.7)
sd(auc.gbh.pers.loop0.75)
sd(auc.gbh.pers.loop0.8)
sd(auc.gbh.pers.loop0.85)
sd(auc.gbh.pers.loop0.9)
sd(auc.gbh.pers.loop0.95)


mean(auc.agp.ab.pers.loop0.5)
mean(auc.agp.ab.pers.loop0.55)
mean(auc.agp.ab.pers.loop0.6)
mean(auc.agp.ab.pers.loop0.65)
mean(auc.agp.ab.pers.loop0.7)
mean(auc.agp.ab.pers.loop0.75)
mean(auc.agp.ab.pers.loop0.8)
mean(auc.agp.ab.pers.loop0.85)
mean(auc.agp.ab.pers.loop0.9)

sd(auc.agp.ab.pers.loop0.55)
sd(auc.agp.ab.pers.loop0.6)
sd(auc.agp.ab.pers.loop0.65)
sd(auc.agp.ab.pers.loop0.7)
sd(auc.agp.ab.pers.loop0.75)
sd(auc.agp.ab.pers.loop0.8)
sd(auc.agp.ab.pers.loop0.85)
sd(auc.agp.ab.pers.loop0.9)


mean(auc.agp.cd.pers.loop0.5)
mean(auc.agp.cd.pers.loop0.55)
mean(auc.agp.cd.pers.loop0.6)
mean(auc.agp.cd.pers.loop0.65)
mean(auc.agp.cd.pers.loop0.7)
mean(auc.agp.cd.pers.loop0.75)
mean(auc.agp.cd.pers.loop0.8)
mean(auc.agp.cd.pers.loop0.85)
mean(auc.agp.cd.pers.loop0.9)

sd(auc.agp.cd.pers.loop0.5)
sd(auc.agp.cd.pers.loop0.55)
sd(auc.agp.cd.pers.loop0.6)
sd(auc.agp.cd.pers.loop0.65)
sd(auc.agp.cd.pers.loop0.7)
sd(auc.agp.cd.pers.loop0.75)
sd(auc.agp.cd.pers.loop0.8)
sd(auc.agp.cd.pers.loop0.85)
sd(auc.agp.cd.pers.loop0.9)

## Make plot

par(mfrow = c(1, 3))

# AGP antibiotics
# Create plotting vectors
groups.agp.ab <- c("0.50", "0.55", "0.60", "0.65", "0.70", "0.75", "0.80", "0.85", "0.90")
means.agp.ab <- c( mean(auc.agp.ab.pers.loop0.5), mean(auc.agp.ab.pers.loop0.55), mean(auc.agp.ab.pers.loop0.6), mean(auc.agp.ab.pers.loop0.65), mean(auc.agp.ab.pers.loop0.7), mean(auc.agp.ab.pers.loop0.75), mean(auc.agp.ab.pers.loop0.8), mean(auc.agp.ab.pers.loop0.85), mean(auc.agp.ab.pers.loop0.9))
sds.agp.ab <-  c( sd(auc.agp.ab.pers.loop0.5), sd(auc.agp.ab.pers.loop0.55), sd(auc.agp.ab.pers.loop0.6), sd(auc.agp.ab.pers.loop0.65), sd(auc.agp.ab.pers.loop0.7), sd(auc.agp.ab.pers.loop0.75), sd(auc.agp.ab.pers.loop0.8), sd(auc.agp.ab.pers.loop0.85), sd(auc.agp.ab.pers.loop0.9))
x_coords.agp.ab <- c(.5, .55, .6, .65, .7, .75, .8, .85, .9) # X-axis positions

# Create other mean and sd vectors to use for setting up plot limits
# AGP C. diff
# groups.agp.cd <- c( "0.60", "0.65", "0.70", "0.75", "0.80", "0.85")
# means.agp.cd <- c( mean(auc.agp.cd.pers.loop0.6), mean(auc.agp.cd.pers.loop0.65), mean(auc.agp.cd.pers.loop0.7), mean(auc.agp.cd.pers.loop0.75), mean(auc.agp.cd.pers.loop0.8), mean(auc.agp.cd.pers.loop0.85))
# sds.agp.cd <-  c( sd(auc.agp.cd.pers.loop0.6), sd(auc.agp.cd.pers.loop0.65), sd(auc.agp.cd.pers.loop0.7), sd(auc.agp.cd.pers.loop0.75), sd(auc.agp.cd.pers.loop0.8), sd(auc.agp.cd.pers.loop0.85))
# x_coords.agp.cd <- c( .6, .65, .7, .75, .8, .85) # X-axis positions

groups.agp.cd <- c("0.50", "0.55", "0.60", "0.65", "0.70", "0.75", "0.80", "0.85", "0.90")
means.agp.cd <- c( mean(auc.agp.cd.pers.loop0.5), mean(auc.agp.cd.pers.loop0.55), mean(auc.agp.cd.pers.loop0.6), mean(auc.agp.cd.pers.loop0.65), mean(auc.agp.cd.pers.loop0.7), mean(auc.agp.cd.pers.loop0.75), mean(auc.agp.cd.pers.loop0.8), mean(auc.agp.cd.pers.loop0.85), mean(auc.agp.cd.pers.loop0.9))
sds.agp.cd <-  c( sd(auc.agp.cd.pers.loop0.5), sd(auc.agp.cd.pers.loop0.55), sd(auc.agp.cd.pers.loop0.6), sd(auc.agp.cd.pers.loop0.65), sd(auc.agp.cd.pers.loop0.7), sd(auc.agp.cd.pers.loop0.75), sd(auc.agp.cd.pers.loop0.8), sd(auc.agp.cd.pers.loop0.85), sd(auc.agp.cd.pers.loop0.9))
x_coords.agp.cd <- c(.5, .55, .6, .65, .7, .75, .8, .85, .9) # X-axis positions

# GBH
groups.gbh <- c("0.60", "0.65", "0.70", "0.75", "0.80", "0.85", "0.90", "0.95")
means.gbh <- c(mean(auc.gbh.pers.loop0.6), mean(auc.gbh.pers.loop0.65), mean(auc.gbh.pers.loop0.7), mean(auc.gbh.pers.loop0.75), mean(auc.gbh.pers.loop0.8), mean(auc.gbh.pers.loop0.85), mean(auc.gbh.pers.loop0.9), mean(auc.gbh.pers.loop0.95))
sds.gbh <-  c(sd(auc.gbh.pers.loop0.6), sd(auc.gbh.pers.loop0.65), sd(auc.gbh.pers.loop0.7), sd(auc.gbh.pers.loop0.75), sd(auc.gbh.pers.loop0.8), sd(auc.gbh.pers.loop0.85), sd(auc.gbh.pers.loop0.9), sd(auc.gbh.pers.loop0.95))
x_coords.gbh <- c(.6, .65, .7 ,.75, .8, .85, .9, .95) # X-axis positions

# Setup plotting area (calculate y limits based on mean +/- sd)
plot(x_coords.agp.ab, means.agp.ab,
     ylim = range(c(means.gbh - sds.gbh, means.gbh + sds.gbh, means.agp.cd - sds.agp.cd, means.agp.cd + sds.agp.cd, means.agp.ab - sds.agp.ab, means.agp.ab + sds.agp.ab)), # Ensure error bars fit
     pch = 19, xlab = "Reference Persistence Threshold", ylab = "Mean AUC +/- 1 SD",
     xaxt = "n", # Turn off default x-axis
     main = "AGP Antibiotics",
     col = "royalblue4", cex = 1.5,
     #log = "x",
     cex.lab = 1.4,
     cex.main = 1.4)

# Add x-axis labels
axis(1, at = x_coords.agp.ab, labels = groups.agp.ab)

# Add error bars using arrows()
# Code 3 indicates arrowheads at both ends; angle 90 makes them flat
arrows(x_coords.agp.ab, means.agp.ab - sds.agp.ab, x_coords.agp.ab, means.agp.ab + sds.agp.ab,
       length = 0.05, angle = 90, code = 3, col = "royalblue4")

# AGP C. diff plot

# Setup plotting area (calculate y limits based on mean +/- sd)
plot(x_coords.agp.cd, means.agp.cd,
     ylim = range(c(means.gbh - sds.gbh, means.gbh + sds.gbh, means.agp.cd - sds.agp.cd, means.agp.cd + sds.agp.cd, means.agp.ab - sds.agp.ab, means.agp.ab + sds.agp.ab)), # Ensure error bars fit
     pch = 19, xlab = "Reference Persistence Threshold", ylab = "Mean AUC +/- 1 SD",
     xaxt = "n", # Turn off default x-axis
     main = "AGP C.diff",
     col = "red4", cex = 1.5,
     #log = "x",
     cex.lab = 1.4,
     cex.main = 1.4)

# Add x-axis labels
axis(1, at = x_coords.agp.cd, labels = groups.agp.cd)

# Add error bars using arrows()
# Code 3 indicates arrowheads at both ends; angle 90 makes them flat
arrows(x_coords.agp.cd, means.agp.cd - sds.agp.cd, x_coords.agp.cd, means.agp.cd + sds.agp.cd,
       length = 0.05, angle = 90, code = 3, col = "red4")

# GBH

# Setup plotting area (calculate y limits based on mean +/- sd)
plot(x_coords.gbh, means.gbh,
     ylim = range(c(means.gbh - sds.gbh, means.gbh + sds.gbh, means.agp.cd - sds.agp.cd, means.agp.cd + sds.agp.cd, means.agp.ab - sds.agp.ab, means.agp.ab + sds.agp.ab)), # Ensure error bars fit
     pch = 19, xlab = "Reference Persistence Threshold", ylab = "Mean AUC +/- 1 SD",
     xaxt = "n", # Turn off default x-axis
     main = "Ghana Breast Health",
     col = "darkgoldenrod4",
     cex = 1.5,
     #log = "x",
     cex.lab = 1.4,
     cex.main = 1.4)

# Add x-axis labels
axis(1, at = x_coords.gbh, labels = groups.gbh)

# Add error bars using arrows()
# Code 3 indicates arrowheads at both ends; angle 90 makes them flat
arrows(x_coords.gbh, means.gbh - sds.gbh, x_coords.gbh, means.gbh + sds.gbh,
       length = 0.05, angle = 90, code = 3, col = "darkgoldenrod4")

