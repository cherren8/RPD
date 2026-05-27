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

# Remove OTU names to leave just numeric matrix
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
healthy.seed <- 1111
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


#### Look at how reference size affects AUC

## Set seed to make analyses reproducible
set.seed(8888)

## AGP dataset:
# Must make ref size smaller here to allow for less overlap between random groups

# Set number of reps
reps.agp <- 200

# set options for test set control size
test.control <- 200
# Set options for control size
control.n.vec <- c(20, 50, 100, 200, 400)

for(control.index in 1:length(control.n.vec)){

  # Create a vector to hold AUCs
  auc.agp.loop <- vector()
  # Create a vector to hold R2 values of explanatory linear model
  r2.agp.loop <- vector()

  for(rep in 1:reps.agp){

    # Set control size
    control.n <- control.n.vec[control.index]

    # Create random control indices
    ref.indices.agp <- sample(1:length(healthy), size = control.n, replace = F)

    # Create reference IDs
    healthy.ref.agp.loop <- healthy[ref.indices.agp]

    # Pull out another set of healthy, plus antibiotic individuals
    ab.test.indices.loop <- c(sample(healthy[-ref.indices.agp], test.control),
                              which(meta.sub1$antibiotic_history == "Week"),
                              which(meta.sub1$antibiotic_history == "Month"))

    # Create reference dataset
    ab.ref.loop <- as.matrix(d.sub1[healthy.ref.agp.loop, ])
    # Create test dataset (contains both controls and cases)
    ab.test.loop <- as.matrix(d.sub1[ab.test.indices.loop, ])
    # Create labels for test dataset
    meta.ab.test.loop <- meta.sub1$antibiotic_history[ab.test.indices.loop]

    # Calculate rpd for antibiotic test dataset
    ab.test.rpd.loop <- calc.rpd(ab.ref.loop, ab.test.loop, ref.pers = 0.7)

    # Create labels to use for GLM
    case.control.binary.abx.loop <- as.numeric(meta.ab.test.loop == "Week" | meta.ab.test.loop == "Month" )

    # Fit the GLM model
    ab.data.loop <- as.data.frame(cbind(case.control.binary.abx.loop, ab.test.rpd.loop))
    model.ab.rpd.loop <- glm(case.control.binary.abx.loop ~ ab.test.rpd.loop , data = ab.data.loop, family = "binomial")

    #summary(model.ab.rpd.loop)

    # Generate predicted probabilities
    preds.ab.rpd.loop <- predict(model.ab.rpd.loop, newdata = ab.data.loop, type = "response")

    # Create the ROC curve and extract AUC
    roc_obj.ab.rpd.loop <- roc(ab.data.loop$case.control.binary.abx.loop, preds.ab.rpd.loop)

    abx.auc.loop <- auc(roc_obj.ab.rpd.loop)

    auc.agp.loop[rep] <- abx.auc.loop

    results.vec.name <- paste0("auc.agp.loop", control.n)

    assign(results.vec.name, auc.agp.loop)

  }

}

mean(auc.agp.loop20)
mean(auc.agp.loop50)
mean(auc.agp.loop100)
mean(auc.agp.loop200)
mean(auc.agp.loop400)

sd(auc.agp.loop20)
sd(auc.agp.loop50)
sd(auc.agp.loop100)
sd(auc.agp.loop200)
sd(auc.agp.loop400)


## Repeat for C diff:

## Set seed to make analyses reproducible
set.seed(8888)

# Set number of reps
reps.agp.cd <- 200
# set options for test set control size
test.control <- 200
# Set options for control size
control.n.vec <- c(20, 50, 100, 200, 400)

for(control.index in 1:length(control.n.vec)){

  # Create a vector to hold AUCs
  auc.agp.cd.loop <- vector()
  # Create a vector to hold R2 values of explanatory linear model
  r2.agp.cd.loop <- vector()

  for(rep in 1:reps.agp.cd){

    # Set control size
    control.n <- control.n.vec[control.index]

    # Create random control indices
    ref.indices.agp.cd <- sample(1:length(healthy), size = control.n, replace = F)

    # Create reference IDs
    healthy.ref.agp.cd.loop <- healthy[ref.indices.agp.cd]

    # Pull out another sample of healthy, plus antibiotic individuals
    cd.test.indices.loop <- c(sample(healthy[-ref.indices.agp.cd], test.control),
                              which(meta.sub1$cdiff == "Diagnosed by a medical professional (doctor, physician assistant)"))

    # Create reference dataset
    cd.ref.loop <- as.matrix(d.sub1[healthy.ref.agp.cd.loop, ])
    # Create test dataset (contains both controls and cases)
    cd.test.loop <- as.matrix(d.sub1[cd.test.indices.loop, ])
    # Create labels for test dataset
    meta.cd.test.loop <- meta.sub1$cdiff[cd.test.indices.loop]

    # Calculate rpd for antibiotic test dataset
    cd.test.rpd.loop <- calc.rpd(cd.ref.loop, cd.test.loop, ref.pers = 0.6)

    # Create labels to use for GLM
    case.control.binary.cd.loop <-  as.numeric(grepl("Diagnosed", meta.cd.test.loop ))

    # Fit the GLM model
    cd.data.loop <- as.data.frame(cbind(case.control.binary.cd.loop, cd.test.rpd.loop))
    model.cd.rpd.loop <- glm(case.control.binary.cd.loop ~ cd.test.rpd.loop , data = cd.data.loop, family = "binomial")

    #summary(model.ab.rpd.loop)

    # Generate predicted probabilities
    preds.cd.rpd.loop <- predict(model.cd.rpd.loop, newdata = cd.data.loop, type = "response")

    # Create the ROC curve and extract AUC
    roc_obj.cd.rpd.loop <- roc(cd.data.loop$case.control.binary.cd.loop, preds.cd.rpd.loop)

    cd.auc.loop <- auc(roc_obj.cd.rpd.loop)

    auc.agp.cd.loop[rep] <- cd.auc.loop

    results.vec.name.cd <- paste0("auc.agp.cd.loop", control.n)

    assign(results.vec.name.cd, auc.agp.cd.loop)

  }

}

mean(auc.agp.cd.loop20)
mean(auc.agp.cd.loop50)
mean(auc.agp.cd.loop100)
mean(auc.agp.cd.loop200)
mean(auc.agp.cd.loop400)

sd(auc.agp.cd.loop20)
sd(auc.agp.cd.loop50)
sd(auc.agp.cd.loop100)
sd(auc.agp.cd.loop200)
sd(auc.agp.cd.loop400)


####################################################################################################
####################################################################################################

## GBH dataset:

# Set number of reps
reps.gbh <- 100

# set options for test set control size
test.control.gbh <- 200
# Set options for control size
control.n.vec.gbh <- c(20, 50, 100, 200) # cannot do 400 because too few samples

set.seed(8888)

for(control.index in 1:length(control.n.vec.gbh)){

  # Create vector to hold AUCs
  auc.ref.size.gbh <- vector()

  for(rep in 1:reps.gbh){

    # Set control size
    control.n <- control.n.vec.gbh[control.index]

    # Pull out all controls
    healthy.gbh <- which(m.sub$case_control__rtl_case_or_contro == "control")

    # Create random control indices
    ref.indices.gbh <- sample(1:length(healthy.gbh), size = control.n, replace = F)

    # Create reference set
    ref <- healthy.gbh[ref.indices.gbh]
    # Create control dataset from sample of remainder of healthy dataset
    control <- sample(healthy.gbh[-ref.indices.gbh], test.control.gbh)
    # Create case dataset from samples labeled "case"
    case <- which(m.sub$case_control__rtl_case_or_contro == "case")

    # Calculate rpd
    rpd.vec.gbh <- calc.rpd(ref = as.matrix(b.num[ref, ]),
                            test = as.matrix(b.num[c(case, control), ]),
                            ref.pers = .7)

    case.control.vec.gbh <- c(rep("case", times = length(case)), rep("control", times = length(control)))

    # Convert case/control vector to 1/0
    case.control.binary.gbh <- as.numeric(case.control.vec.gbh == "case")

    gbh.data.ref.pers <- as.data.frame(cbind(case.control.binary.gbh, rpd.vec.gbh))
    model.br.rpd <- glm(case.control.binary.gbh ~ rpd.vec.gbh , data = gbh.data.ref.pers, family = "binomial")

    # Generate predicted probabilities
    preds.br.rpd <- predict(model.br.rpd, newdata = gbh.data.ref.pers, type = "response")

    # Create and plot the ROC curve
    roc_obj.br.rpd <- roc(gbh.data.ref.pers$case.control.binary.gbh, preds.br.rpd)

    auc.ref.size.gbh[rep] <- auc(roc_obj.br.rpd)

  }

  results.size.name.gbh <- paste0("auc.gbh.size.loop", control.n)

  assign(results.size.name.gbh, auc.ref.size.gbh)

  hist(auc.ref.size.gbh)

}

# Print all Mean and SDs
mean(auc.agp.loop20)
mean(auc.agp.loop50)
mean(auc.agp.loop100)
mean(auc.agp.loop200)
mean(auc.agp.loop400)

sd(auc.agp.loop20)
sd(auc.agp.loop50)
sd(auc.agp.loop100)
sd(auc.agp.loop200)
sd(auc.agp.loop400)

mean(auc.agp.cd.loop20)
mean(auc.agp.cd.loop50)
mean(auc.agp.cd.loop100)
mean(auc.agp.cd.loop200)
mean(auc.agp.cd.loop400)

sd(auc.agp.cd.loop20)
sd(auc.agp.cd.loop50)
sd(auc.agp.cd.loop100)
sd(auc.agp.cd.loop200)
sd(auc.agp.cd.loop400)

mean(auc.gbh.size.loop20)
mean(auc.gbh.size.loop50)
mean(auc.gbh.size.loop100)
mean(auc.gbh.size.loop200)

sd(auc.gbh.size.loop20)
sd(auc.gbh.size.loop50)
sd(auc.gbh.size.loop100)
sd(auc.gbh.size.loop200)


# Make plotting values so plot area can be set up
# AGP Antibiotics
par(mfrow = c(1, 3))
#  Create plotting vectors
groups.agp.ab <- c("20", "50", "100", "200", "400")
means.agp.ab <- c(mean(auc.agp.loop20), mean(auc.agp.loop50), mean(auc.agp.loop100), mean(auc.agp.loop200), mean(auc.agp.loop400))
sds.agp.ab <- c(sd(auc.agp.loop20), sd(auc.agp.loop50), sd(auc.agp.loop100), sd(auc.agp.loop200), sd(auc.agp.loop400))
x_coords.agp.ab <- c(20, 50, 100, 200, 400) # X-axis positions

#  AGP C. diff
groups.agp.cd <- c("20", "50", "100", "200", "400")
means.agp.cd <- c(mean(auc.agp.cd.loop20), mean(auc.agp.cd.loop50), mean(auc.agp.cd.loop100), mean(auc.agp.cd.loop200), mean(auc.agp.cd.loop400))
sds.agp.cd <- c(sd(auc.agp.cd.loop20), sd(auc.agp.cd.loop50), sd(auc.agp.cd.loop100), sd(auc.agp.cd.loop200), sd(auc.agp.cd.loop400))
x_coords.agp.cd <- c(20, 50, 100, 200, 400) # X-axis positions

#  GBH
groups.gbh <- c("20", "50", "100", "200")
means.gbh <- c(mean(auc.gbh.size.loop20), mean(auc.gbh.size.loop50), mean(auc.gbh.size.loop100), mean(auc.gbh.size.loop200))
sds.gbh <- c(sd(auc.gbh.size.loop20), sd(auc.gbh.size.loop50), sd(auc.gbh.size.loop100), sd(auc.gbh.size.loop200))
x_coords.gbh <- c(20, 50, 100, 200) # X-axis positions

#  Plot data
plot(x_coords.agp.ab, means.agp.ab,
     ylim = range(c(means.agp.ab - sds.agp.ab, means.agp.ab + sds.agp.ab, means.agp.cd - sds.agp.cd, means.agp.cd + sds.agp.cd)), # Ensure error bars fit
     pch = 19, xlab = "Reference Dataset Samples", ylab = "Mean AUC +/- 1 SD",
     xaxt = "n", # Turn off default x-axis
     main = "AGP Antibiotics",
     col = "royalblue4", cex = 1.5,
     log = "x", xlim = c(15, 500),
     cex.lab = 1.4,
     cex.main = 1.4)

#  Add x-axis labels
axis(1, at = x_coords.agp.ab, labels = groups.agp.ab)

#  Add error bars using arrows()
# Code 3 indicates arrowheads at both ends; angle 90 makes them flat
arrows(x_coords.agp.ab, means.agp.ab - sds.agp.ab, x_coords.agp.ab, means.agp.ab + sds.agp.ab,
       length = 0.05, angle = 90, code = 3, col = "royalblue4")

# AGP C. diff

#  Plot data
plot(x_coords.agp.cd, means.agp.cd,
     ylim = range(c(means.agp.ab - sds.agp.ab, means.agp.ab + sds.agp.ab, means.agp.cd - sds.agp.cd, means.agp.cd + sds.agp.cd)), # Ensure error bars fit
     pch = 19, xlab = "Reference Dataset Samples", ylab = "Mean AUC +/- 1 SD",
     xaxt = "n", # Turn off default x-axis
     main = "AGP C. Diff",
     col = "red4", cex = 1.5,
     log = "x", xlim = c(15, 500),
     cex.lab = 1.4,
     cex.main = 1.4)

#  Add x-axis labels
axis(1, at = x_coords.agp.cd, labels = groups.agp.cd)

#  Add error bars using arrows()
# Code 3 indicates arrowheads at both ends; angle 90 makes them flat
arrows(x_coords.agp.cd, means.agp.cd - sds.agp.cd, x_coords.agp.cd, means.agp.cd + sds.agp.cd,
       length = 0.05, angle = 90, code = 3, col = "red4")


# GBH

#  Plot data
plot(x_coords.gbh, means.gbh,
     ylim = range(c(means.agp.ab - sds.agp.ab, means.agp.ab + sds.agp.ab, means.agp.cd - sds.agp.cd, means.agp.cd + sds.agp.cd)), # Ensure error bars fit
     pch = 19, xlab = "Reference Dataset Samples", ylab = "Mean AUC +/- 1 SD",
     xaxt = "n", # Turn off default x-axis
     main = "Ghana Breast Health",
     col = "darkgoldenrod4", cex = 1.5,
     log = "x", xlim = c(15, 250),
     cex.lab = 1.4,
     cex.main = 1.4)

#  Add x-axis labels
axis(1, at = x_coords.gbh, labels = groups.gbh)

#  Add error bars using arrows()
# Code 3 indicates arrowheads at both ends; angle 90 makes them flat
arrows(x_coords.gbh, means.gbh - sds.gbh, x_coords.gbh, means.gbh + sds.gbh,
       length = 0.05, angle = 90, code = 3, col = "darkgoldenrod4")

