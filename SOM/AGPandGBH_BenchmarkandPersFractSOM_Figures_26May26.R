#### Supplementary Figures for RPD Methods Paper. Begins with script used in main text:

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
bc.ref.mean.gbh <- apply(b.num[ref,], 2, mean)
# Create matrix of case and control samples for use in Bray-Curtis analysis
b.num.bc <- rbind(b.num[case, ], b.num[control, ])

# Calculate BCD to centroid for each sample
bc.vec.ref <- as.matrix(vegdist(rbind(bc.ref.mean.gbh, b.num.bc)))[-1, 1]

# Implement same logistic regression
summary(glm(as.factor(case.control.binary.gbh) ~ bc.vec.ref, family = "binomial"))

## Calculate diversity metrics
ref.div.gbh <- vegan::diversity(b.num[ref, ])
control.div.gbh <- vegan::diversity(b.num[control, ])
case.div.gbh <- vegan::diversity(b.num[case, ])
case.control.div.gbh <- c(case.div.gbh, control.div.gbh)

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
ab.ref.mean <- apply(ab.ref, 2, mean)
# Calculate BCD to centroid
bc.vec.ab <- as.matrix(vegdist(rbind(ab.ref.mean, ab.test)))[-1, 1]

# Repeat logistic regression
summary(glm(as.factor(case.control.binary.abx) ~ bc.vec.ab, family = "binomial"))

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
cd.ref.mean <- apply(cd.ref, 2, mean)

bc.vec.cd <- as.matrix(vegdist(rbind(cd.ref.mean, cd.test)))[-1, 1]

# Implement logistic regression
summary(glm(as.factor(case.control.binary.cd) ~ bc.vec.cd, family = "binomial"))

####################################################################################################
####################################################################################################

# Make figures of rpd versus diversity and BCD

# First, do a check to make sure that sequencing depth does not imapct diversity too much:
# Calculate raw diversity, then remove all values below the threshold of the smallest value consistently observed in the data, and re-calculate
# If correlation of diversity metrics in the two scenarios is high, then sequencing depth does not affect diversity

d.div <- vegan::diversity(d.sub1)
#
# # Plot min detection in each sample
# min.non.zero <- function(vec){
#   min.pos <- min(vec[vec>0])
# }
#
# # Find detection threshold by looking at the smallest observed values in each sample
# hist(log(apply(d, 1, min.non.zero)))
# hist(log(apply(d[cd.test.indices,], 1, min.non.zero)))
# exp(-8) # this value is pretty consistently detected in all datsets
#
# # Create a dataset removing all values smaller than e^-8
# d.cut8 <- d.sub1
# d.cut8[d.cut8 < exp(-8)] <- 0
#
# # Calculate diversity on this thresholded dataset
# d8.div <- diversity(d.cut8)
#
# # Look at correlation of diversity calculated on two datasets
# plot(d.div, d8.div) # excellent, cutting out small values has almost no impact
# cor(d.div, d8.div) # 0.998. Valid to use raw diversity

####
#### Make rpd vs diversity and rpd vs BCD figures:

# par(mfrow = c(3, 2))
# plot(ab.test.rpd ~ d.div[ab.test.indices],
#      ylim = range(c(ab.test.rpd, cd.test.rpd, rpd.vec.gbh + .02), na.rm =  T),
#      xlim = range(c(d.div, case.control.div.gbh)),
#      col = c(rep(adjustcolor("royalblue1", alpha = .5), times = sum(case.control.binary.abx == 0) ), rep(adjustcolor("royalblue4", alpha = .5), times = sum(case.control.binary.abx == 1) ) ),
#      xlab = "Diversity",
#      ylab = "Sample RPD",
#      main = "Antibiotic History",
#      pch = 16,
#      cex= 1.2)
# text(labels = paste("r = ", round(cor(ab.test.rpd , d.div[ab.test.indices], use = "pairwise.complete.obs"), digits = 3) ) , x = 1, y = .35, cex = 1.2)
#
# legend("bottomleft", legend = c("control", "case"),
#        col = c(adjustcolor("royalblue1", alpha = .5), adjustcolor("royalblue4", alpha = .5)),
#        pch = 16,
#        cex = 1.2)
#
# cor(ab.test.rpd , d.div[ab.test.indices], use = "pairwise.complete.obs")
#
# plot(ab.test.rpd ~ bc.vec.ab,
#      xlim = range(c(bc.vec.ab, bc.vec.cd, bc.vec.ref), na.rm =  T),
#      ylim = range(c(ab.test.rpd, cd.test.rpd, rpd.vec.gbh + .02), na.rm =  T),
#      col = c(rep(adjustcolor("royalblue1", alpha = .5), times = sum(case.control.binary.abx == 0) ), rep(adjustcolor("royalblue4", alpha = .5), times = sum(case.control.binary.abx == 1) ) ),
#      xlab = "BCD to Ref. Centroid",
#      ylab = "Sample RPD",
#      main = "Antibiotic History",
#      pch = 16,
#      cex= 1.2)
# text(labels = paste("r = ", round(cor(ab.test.rpd , bc.vec.ab, use = "pairwise.complete.obs"), digits = 3) ) , x = .5, y = .35, cex = 1.2)
#
# cor(ab.test.rpd , bc.vec.ab, use = "pairwise.complete.obs")
#
# plot(cd.test.rpd ~ d.div[cd.test.indices],
#      ylim = range(c(ab.test.rpd, cd.test.rpd, rpd.vec.gbh + .02), na.rm =  T),
#      xlim = range(c(d.div, case.control.div.gbh)),
#      col = c(rep(adjustcolor("red1", alpha = .4), times = sum(case.control.binary.cd == 0) ), rep(adjustcolor("red4", alpha = .4), times = sum(case.control.binary.cd == 1) ) ),
#      xlab = "Diversity",
#      ylab = "Sample RPD",
#      main = "C. diff Diagnosis",
#      pch = 16,
#      cex= 1.2)
# text(labels = paste("r = ", round(cor(cd.test.rpd , d.div[cd.test.indices], use = "pairwise.complete.obs"), digits = 3) ) , x = 1, y = .35, cex = 1.2)
# legend("bottomleft", legend = c("control", "case"),
#        col = c(adjustcolor("red1", alpha = .4), adjustcolor("red4", alpha = .4)),
#        pch = 16,
#        cex = 1.2)
#
# cor(cd.test.rpd , d.div[cd.test.indices], use = "pairwise.complete.obs")
#
# plot(cd.test.rpd ~ bc.vec.cd,
#      xlim = range(c(bc.vec.ab, bc.vec.cd, bc.vec.ref), na.rm =  T),
#      ylim = range(c(ab.test.rpd, cd.test.rpd, rpd.vec.gbh + .02), na.rm =  T),
#      col = c(rep(adjustcolor("red1", alpha = .4), times = sum(case.control.binary.cd == 0) ), rep(adjustcolor("red4", alpha = .4), times = sum(case.control.binary.cd == 1) ) ),
#      xlab = "BCD to Ref. Centroid",
#      ylab = "Sample RPD",
#      main = "C.diff Diagnosis",
#      pch = 16,
#      cex= 1.2)
# text(labels = paste("r = ", round(cor(cd.test.rpd , bc.vec.cd, use = "pairwise.complete.obs"), digits = 3) ) , x = .5, y = .35, cex = 1.2)
#
# cor(cd.test.rpd , bc.vec.cd, use = "pairwise.complete.obs")
#
# plot(rpd.vec.gbh ~ case.control.div.gbh,
#      ylim = range(c(ab.test.rpd, cd.test.rpd, rpd.vec.gbh + .02), na.rm =  T),
#      xlim = range(c(d.div, case.control.div.gbh)),
#      col = c(rep(adjustcolor("darkgoldenrod4", alpha = .4), times = sum(case.control.vec.gbh == "case") ), rep(adjustcolor("darkgoldenrod1", alpha = .4), times = sum(case.control.vec.gbh == "control") ) ),
#      xlab = "Diversity",
#      ylab = "Sample RPD",
#      main = "Breast Health",
#      pch = 16,
#      cex= 1.2)
# text(labels = paste("r = ", round(cor(rpd.vec.gbh , case.control.div.gbh, use = "pairwise.complete.obs"), digits = 3) ) , x = 1, y = .35, cex = 1.2)
#
# legend("bottomleft", legend = c("control", "case"),
#        col = c(adjustcolor("darkgoldenrod1", alpha = .4), adjustcolor("darkgoldenrod4", alpha = .4)),
#        pch = 16,
#        cex = 1.2)
#
# cor(rpd.vec.gbh , case.control.div.gbh, use = "pairwise.complete.obs")
#
# plot(rpd.vec.gbh ~ bc.vec.ref,
#      xlim = range(c(bc.vec.ab, bc.vec.cd, bc.vec.ref), na.rm =  T),
#      ylim = range(c(ab.test.rpd, cd.test.rpd, rpd.vec.gbh + .02), na.rm =  T),
#      col = c(rep(adjustcolor("darkgoldenrod4", alpha = .4), times = sum(case.control.vec.gbh == "case") ), rep(adjustcolor("darkgoldenrod1", alpha = .4), times = sum(case.control.vec.gbh == "control") ) ),
#      xlab = "BCD to Ref. Centroid",
#      ylab = "Sample RPD",
#      main = "Breast Health",
#      pch = 16,
#      cex= 1.2)
# text(labels = paste("r = ", round(cor(rpd.vec.gbh , bc.vec.ref, use = "pairwise.complete.obs"), digits = 3) ) , x = .5, y = .35, cex = 1.2)
#
# cor(rpd.vec.gbh , bc.vec.ref, use = "pairwise.complete.obs")
#
# ## Now make violin plots
# par(mfrow = c(1, 3))
#
# # Relevel label vectors to make plots in correct order
# meta.cd.test <- relevel(as.factor(meta.cd.test), ref = "I do not have this condition")
#
# case.control.vec.gbh <- relevel(as.factor(case.control.vec.gbh), ref = "control")
#
# par(cex.lab = 1.3)
#
# vioplot(ab.test.rpd ~ case.control.binary.abx,
#         names = c("None", "Recent Antibiotics"),
#         xlab = "Antibiotic Status",
#         ylab = "Sample RPD",
#         col = c(adjustcolor("royalblue1", alpha = .5), adjustcolor("royalblue4", alpha = .5)) ,
#         cex.axis = 1.2 )
#
# vioplot(cd.test.rpd ~ meta.cd.test,
#         names = c( "No diagnosis", "Positive C.diff diagnosis"),
#         xlab = "C. diff diagnosis Status",
#         ylab = "Sample RPD",
#         cex.axis = 1.2 ,
#         col = c(adjustcolor("red1", alpha = .5), adjustcolor("red4", alpha = .5)) )
#
# vioplot(rpd.vec.gbh ~ case.control.vec.gbh,
#         names = c("Control", "Disease"),
#         xlab = "Breast Disease Status",
#         ylab = "Sample RPD",
#         cex.axis = 1.2 ,
#         col = c(adjustcolor("darkgoldenrod1", alpha = .5), adjustcolor("darkgoldenrod4", alpha = .5)) )


######################################################################################
######################################################################################
#### Conduct additional benchmarks and add to ROC curves

#### AGP datasets Rank Shift Analysis

# Find ranks of centroid, keeping only taxa that meet persistence threshold
rank.mean.vec.agp <- rank(ab.ref.mean[colMeans(ceiling(ab.ref)) > .7 ])
length(rank.mean.vec.agp)

# Make new dataset with only the taxa present at the persistence threshold
ab.test.sub <- ab.test[, colMeans(ceiling(ab.ref)) > .7]
dim(ab.test.sub) # should be same number of columns as length of rank.mean.vec

# Make new dataset with only the taxa present at the persistence threshold
cd.test.sub <- cd.test[, colMeans(ceiling(ab.ref)) > .7]
dim(cd.test.sub)

# Rank the taxa in each sample (each row)
ref.mat.rank.ab <- apply(ab.test.sub, 1, rank)

ref.mat.rank.cd <- apply(cd.test.sub, 1, rank)

# Subtract the ranks in the centroid from the ranks in the samples
rank.diffs.ab <- sweep(ref.mat.rank.ab, 1, rank.mean.vec.agp)

rank.diffs.cd <- sweep(ref.mat.rank.cd, 1, rank.mean.vec.agp)

# Find the mean absolute rank difference
rank.diff.means.ab <- apply(abs(rank.diffs.ab), 2, mean)

rank.diff.means.cd <- apply(abs(rank.diffs.cd), 2, mean)

#### AGP datasets Aitchison distance with robust clr. Use both subset and full dataset
ait.fract <- .7

## AGP Ab

ait.dist.ab.ref <- as.matrix(vegdist(rbind(ab.ref.mean[colMeans(ceiling(ab.ref)) > ait.fract], ab.test[, colMeans(ceiling(ab.ref)) > ait.fract]), method = "robust.aitchison"))[-1, 1]

ait.dist.ab.ref.full <- as.matrix(vegdist(rbind(ab.ref.mean, ab.test), method = "robust.aitchison"))[-1, 1]


## AGP C. diff

ait.dist.cd.ref <- as.matrix(vegdist(rbind(cd.ref.mean[colMeans(ceiling(cd.ref)) > ait.fract], cd.test[, colMeans(ceiling(ab.ref)) > ait.fract]), method = "robust.aitchison"))[-1, 1]

ait.dist.cd.ref.full <- as.matrix(vegdist(rbind(cd.ref.mean, cd.test), method = "robust.aitchison"))[-1, 1]

#### Euclidean distance on subset

euc.dist.ab.ref <- as.matrix(vegdist(rbind(ab.ref.mean[colMeans(ceiling(ab.ref)) > .7], ab.test[, colMeans(ceiling(ab.ref)) > .7]), method = "euclidean"))[-1, 1]

euc.dist.cd.ref <- as.matrix(vegdist(rbind(cd.ref.mean[colMeans(ceiling(cd.ref)) > .7], cd.test[, colMeans(ceiling(ab.ref)) > .7]), method = "euclidean"))[-1, 1]

#### Median Bray-Curtis to references
# This manual way produces the same as subsetting the full matrix, just 100x slower
# med.bc.ab <- vector()
# for(ab.sample in 1:dim(ab.test)[1]){
#   med.bc.ab[ab.sample] <- median(as.matrix(vegdist(rbind(ab.test[ab.sample,], ab.ref)))[-1, 1])
# }
#
# med.bc.cd <- vector()
# for(cd.sample in 1:dim(cd.test)[1]){
#   med.bc.cd[cd.sample] <- median(as.matrix(vegdist(rbind(cd.test[cd.sample,], cd.ref)))[-1, 1])
# }

full.ab.bc.mat <- as.matrix(vegdist(rbind(ab.test, ab.ref)))
full.ab.bc.mat[1:4, 1:4]
ab.bc.toref.mat <- full.ab.bc.mat[(1:dim(ab.test)[1]), -c(1:dim(ab.test)[1])]
dim(ab.bc.toref.mat)
ab.bc.toref.mat[1:4, 1:4]
med.bc.ab <- apply(ab.bc.toref.mat, 1, median)

full.cd.bc.mat <- as.matrix(vegdist(rbind(cd.test, cd.ref)))
full.cd.bc.mat[1:4, 1:4]
cd.bc.toref.mat <- full.cd.bc.mat[(1:dim(cd.test)[1]), -c(1:dim(cd.test)[1])]
dim(cd.bc.toref.mat)
cd.bc.toref.mat[1:4, 1:4]
med.bc.cd <- apply(cd.bc.toref.mat, 1, median)

#### BCD on subset

bc.vec.ab.sub <- as.matrix(vegdist(rbind(ab.ref.mean[colMeans(ceiling(ab.ref)) > .7], ab.test[, colMeans(ceiling(ab.ref)) > .7])))[-1, 1]

bc.vec.cd.sub <- as.matrix(vegdist(rbind(cd.ref.mean[colMeans(ceiling(cd.ref)) > .7], cd.test[, colMeans(ceiling(cd.ref)) > .7])))[-1, 1]


#### ROC Curves
par(mfrow = c(1, 1))

par(cex.lab = 1.3, cex.axis = 1.1)

## ABX usage

# 1. Fit the GLM model
ab.data <- as.data.frame(cbind(case.control.binary.abx, ab.test.rpd, bc.vec.ab, rank.diff.means.ab, ait.dist.ab.ref, ait.dist.ab.ref.full, euc.dist.ab.ref, med.bc.ab, bc.vec.ab.sub))

model.ab.rpd <- glm(case.control.binary.abx ~ ab.test.rpd , data = ab.data, family = "binomial")
model.ab.bcd <- glm(case.control.binary.abx ~ bc.vec.ab , data = ab.data, family = "binomial")
model.ab.rank <- glm(case.control.binary.abx ~ rank.diff.means.ab , data = ab.data, family = "binomial")
model.ab.ait <- glm(case.control.binary.abx ~ ait.dist.ab.ref , data = ab.data, family = "binomial")
model.ab.aitf <- glm(case.control.binary.abx ~ ait.dist.ab.ref.full , data = ab.data, family = "binomial")
model.ab.euc <- glm(case.control.binary.abx ~ euc.dist.ab.ref , data = ab.data, family = "binomial")
model.ab.bcdm <- glm(case.control.binary.abx ~ med.bc.ab , data = ab.data, family = "binomial")
model.ab.bcds <- glm(case.control.binary.abx ~ bc.vec.ab.sub , data = ab.data, family = "binomial")

# 2. Generate predicted probabilities (type = "response" is crucial)
preds.ab.rpd <- predict(model.ab.rpd, newdata = ab.data, type = "response")
preds.ab.bcd <- predict(model.ab.bcd, newdata = ab.data, type = "response")
preds.ab.rank <- predict(model.ab.rank, newdata = ab.data, type = "response")
preds.ab.ait <- predict(model.ab.ait, newdata = ab.data, type = "response")
preds.ab.aitf <- predict(model.ab.aitf, newdata = ab.data, type = "response")
preds.ab.euc <- predict(model.ab.euc, newdata = ab.data, type = "response")
preds.ab.bcdm <- predict(model.ab.bcdm, newdata = ab.data, type = "response")
preds.ab.bcds <- predict(model.ab.bcds, newdata = ab.data, type = "response")

# 3. Create and plot the ROC curve
roc_obj.ab.rpd <- roc(ab.data$case.control.binary.abx, preds.ab.rpd)
plot(roc_obj.ab.rpd, print.auc = TRUE, col = "blue", main = "ROC for AGP Antibiotics",
     print.auc.x = .4,
     print.auc.y = .55,
     print.auc.cex = 1.2)

roc_obj.ab.bcd <- roc(ab.data$case.control.binary, preds.ab.bcd)
plot(roc_obj.ab.bcd, print.auc = TRUE, col = "red",
     add = T,
     print.auc.x = .4,
     print.auc.y = .47,
     print.auc.cex = 1.2)

roc_obj.ab.bcdm <- roc(ab.data$case.control.binary, preds.ab.bcdm)
plot(roc_obj.ab.bcdm, print.auc = TRUE, col = "maroon",
     add = T,
     print.auc.x = .4,
     print.auc.y = .39,
     print.auc.cex = 1.2)

roc_obj.ab.bcds <- roc(ab.data$case.control.binary, preds.ab.bcds)
plot(roc_obj.ab.bcds, print.auc = TRUE, col = "purple",
     add = T,
     print.auc.x = .4,
     print.auc.y = .31,
     print.auc.cex = 1.2)

roc_obj.ab.ait <- roc(ab.data$case.control.binary, preds.ab.ait)
plot(roc_obj.ab.ait, print.auc = TRUE, col = "yellowgreen",
     add = T,
     print.auc.x = .4,
     print.auc.y = .23,
     print.auc.cex = 1.2)

roc_obj.ab.aitf <- roc(ab.data$case.control.binary, preds.ab.aitf)
plot(roc_obj.ab.aitf, print.auc = TRUE, col = "forestgreen",
     add = T,
     print.auc.x = .4,
     print.auc.y = .15,
     print.auc.cex = 1.2)

roc_obj.ab.euc <- roc(ab.data$case.control.binary, preds.ab.euc)
plot(roc_obj.ab.euc, print.auc = TRUE, col = "orange",
     add = T,
     print.auc.x = .4,
     print.auc.y = .07,
     print.auc.cex = 1.2)

roc_obj.ab.rank <- roc(ab.data$case.control.binary, preds.ab.rank)
plot(roc_obj.ab.rank, print.auc = TRUE, col = "magenta",
     add = T,
     print.auc.x = .4,
     print.auc.y = -0.01,
     print.auc.cex = 1.2)

legend("topleft", legend = c("RPD", "BCD","BCDM", "BCDS", "AIT", "AITF", "EUC", "RS"), col = c("blue", "red", "maroon", "purple", "yellowgreen", "forestgreen", "orange", "magenta"), lwd = 2, bty = "n")

## C. diff
# 1. Fit the GLM model
cd.data <- as.data.frame(cbind(case.control.binary.cd, cd.test.rpd, bc.vec.cd, rank.diff.means.cd, ait.dist.cd.ref, ait.dist.cd.ref.full, euc.dist.cd.ref, med.bc.cd, bc.vec.cd.sub))

model.cd.rpd <- glm(case.control.binary.cd ~ cd.test.rpd , data = cd.data, family = "binomial")
model.cd.bcd <- glm(case.control.binary.cd ~ bc.vec.cd , data = cd.data, family = "binomial")
model.cd.rank <- glm(case.control.binary.cd ~ rank.diff.means.cd , data = cd.data, family = "binomial")
model.cd.ait <- glm(case.control.binary.cd ~ ait.dist.cd.ref , data = cd.data, family = "binomial")
model.cd.aitf <- glm(case.control.binary.cd ~ ait.dist.cd.ref.full , data = cd.data, family = "binomial")
model.cd.euc <- glm(case.control.binary.cd ~ euc.dist.cd.ref , data = cd.data, family = "binomial")
model.cd.bcdm <- glm(case.control.binary.cd ~ med.bc.cd , data = cd.data, family = "binomial")
model.cd.bcds <- glm(case.control.binary.cd ~ bc.vec.cd.sub , data = cd.data, family = "binomial")

# 2. Generate predicted probcdilities (type = "response" is crucial)
preds.cd.rpd <- predict(model.cd.rpd, newdata = cd.data, type = "response")
preds.cd.bcd <- predict(model.cd.bcd, newdata = cd.data, type = "response")
preds.cd.rank <- predict(model.cd.rank, newdata = cd.data, type = "response")
preds.cd.ait <- predict(model.cd.ait, newdata = cd.data, type = "response")
preds.cd.aitf <- predict(model.cd.aitf, newdata = cd.data, type = "response")
preds.cd.euc <- predict(model.cd.euc, newdata = cd.data, type = "response")
preds.cd.bcdm <- predict(model.cd.bcdm, newdata = cd.data, type = "response")
preds.cd.bcds <- predict(model.cd.bcds, newdata = cd.data, type = "response")

# 3. Create and plot the ROC curve
roc_obj.cd.rpd <- roc(cd.data$case.control.binary.cd, preds.cd.rpd)
plot(roc_obj.cd.rpd, print.auc = TRUE, col = "blue", main = "ROC for AGP C. diff",
     print.auc.x = .4,
     print.auc.y = .55,
     print.auc.cex = 1.2)

roc_obj.cd.bcd <- roc(cd.data$case.control.binary, preds.cd.bcd)
plot(roc_obj.cd.bcd, print.auc = TRUE, col = "red",
     add = T,
     print.auc.x = .4,
     print.auc.y = .47,
     print.auc.cex = 1.2)

roc_obj.cd.bcdm <- roc(cd.data$case.control.binary, preds.cd.bcdm)
plot(roc_obj.cd.bcdm, print.auc = TRUE, col = "maroon",
     add = T,
     print.auc.x = .4,
     print.auc.y = .39,
     print.auc.cex = 1.2)

roc_obj.cd.bcds <- roc(cd.data$case.control.binary, preds.cd.bcds)
plot(roc_obj.cd.bcds, print.auc = TRUE, col = "purple",
     add = T,
     print.auc.x = .4,
     print.auc.y = .31,
     print.auc.cex = 1.2)

roc_obj.cd.ait <- roc(cd.data$case.control.binary, preds.cd.ait)
plot(roc_obj.cd.ait, print.auc = TRUE, col = "yellowgreen",
     add = T,
     print.auc.x = .4,
     print.auc.y = .23,
     print.auc.cex = 1.2)

roc_obj.cd.aitf <- roc(cd.data$case.control.binary, preds.cd.aitf)
plot(roc_obj.cd.aitf, print.auc = TRUE, col = "forestgreen",
     add = T,
     print.auc.x = .4,
     print.auc.y = .15,
     print.auc.cex = 1.2)

roc_obj.cd.euc <- roc(cd.data$case.control.binary, preds.cd.euc)
plot(roc_obj.cd.euc, print.auc = TRUE, col = "orange",
     add = T,
     print.auc.x = .4,
     print.auc.y = .07,
     print.auc.cex = 1.2)

roc_obj.cd.rank <- roc(cd.data$case.control.binary, preds.cd.rank)
plot(roc_obj.cd.rank, print.auc = TRUE, col = "magenta",
     add = T,
     print.auc.x = .4,
     print.auc.y = -0.01,
     print.auc.cex = 1.2)

legend("topleft", legend = c("RPD", "BCD","BCDM", "BCDS", "AIT", "AITF", "EUC", "RS"), col = c("blue", "red", "maroon", "purple", "yellowgreen", "forestgreen", "orange", "magenta"), lwd = 2, bty = "n")


## Breast Health

#### Conduct rank shift analysis with same cutoff

# Find ranks of centroid, keeping only taxa that meet persistence threshold
rank.mean.vec.gbh <- rank(bc.ref.mean.gbh[colMeans(ceiling(b.num[ref, ])) > .7 ])
length(rank.mean.vec.gbh)

# Make a new dataset with only the taxa present at the persistence threshold
gbh.test.sub <- b.num.bc[, colMeans(ceiling(b.num[ref, ])) > .7]
dim(gbh.test.sub) # should be same number of columns as length of rank.mean.vec.gbh

# Rank the taxa in each sample (each row)
ref.mat.rank.gbh <- apply(gbh.test.sub, 1, rank)

# Subtract the ranks in the centroid from the ranks in the samples
rank.diffs.gbh <- sweep(ref.mat.rank.gbh, 1, rank.mean.vec.gbh)

# Find the mean absolute rank difference
rank.diff.means.gbh <- apply(abs(rank.diffs.gbh), 2, mean)

#### Aitchison distance on subset and full dataset

ait.dist.gbh.ref <- as.matrix(vegdist(rbind(bc.ref.mean.gbh[colMeans(ceiling(b.num[ref, ])) > .7], b.num.bc[, colMeans(ceiling(b.num[ref, ])) > .7]), method = "robust.aitchison"))[-1, 1]

ait.dist.gbh.ref.full <- as.matrix(vegdist(rbind(bc.ref.mean.gbh, b.num.bc), method = "robust.aitchison"))[-1, 1]

#### Euclidean distance on subset
euc.dist.gbh.ref <- as.matrix(vegdist(rbind(bc.ref.mean.gbh[colMeans(ceiling(b.num[ref, ])) > .7], b.num.bc[, colMeans(ceiling(b.num[ref, ])) > .7]), method = "euclidean"))[-1, 1]

#### Median BCD to reference samples

full.br.bc.mat <- as.matrix(vegdist(rbind(b.num[c(case, control), ], b.num[ref, ])))
dim(full.br.bc.mat)
full.br.bc.mat[1:4, 1:4]
br.bc.toref.mat <- full.br.bc.mat[(1:length(c(case, control))), -c(1:length(c(case, control)))]
dim(br.bc.toref.mat)
cd.bc.toref.mat[1:4, 1:4]
med.bc.br <- apply(br.bc.toref.mat, 1, median)

#### BCD on subset

bc.vec.ref.sub <- as.matrix(vegdist(rbind(bc.ref.mean.gbh[colMeans(ceiling(b.num[ref, ])) > .7], b.num.bc[, colMeans(ceiling(b.num[ref, ])) > .7])))[-1, 1]

# 1. Fit the GLM model
br.data <- as.data.frame(cbind(case.control.binary.gbh, rpd.vec.gbh, bc.vec.ref, rank.diff.means.gbh, ait.dist.gbh.ref, ait.dist.gbh.ref.full, euc.dist.gbh.ref, med.bc.br, bc.vec.ref.sub))
model.br.rpd <- glm(case.control.binary.gbh ~ rpd.vec.gbh , data = br.data, family = "binomial")
model.br.bcd <- glm(case.control.binary.gbh ~ bc.vec.ref , data = br.data, family = "binomial")
model.br.rank <- glm(case.control.binary.gbh ~ rank.diff.means.gbh , data = br.data, family = "binomial")
model.br.ait <- glm(case.control.binary.gbh ~ ait.dist.gbh.ref , data = br.data, family = "binomial")
model.br.aitf <- glm(case.control.binary.gbh ~ ait.dist.gbh.ref.full , data = br.data, family = "binomial")
model.br.euc <- glm(case.control.binary.gbh ~ euc.dist.gbh.ref , data = br.data, family = "binomial")
model.br.bcdm <- glm(case.control.binary.gbh ~ med.bc.br , data = br.data, family = "binomial")
model.br.bcds <- glm(case.control.binary.gbh ~ bc.vec.ref.sub , data = br.data, family = "binomial")

# 2. Generate predicted probabilities (type = "response" is crucial)
preds.br.rpd <- predict(model.br.rpd, newdata = br.data, type = "response")
preds.br.bcd <- predict(model.br.bcd, newdata = br.data, type = "response")
preds.br.rank <- predict(model.br.rank, newdata = br.data, type = "response")
preds.br.ait <- predict(model.br.ait, newdata = br.data, type = "response")
preds.br.aitf <- predict(model.br.aitf, newdata = br.data, type = "response")
preds.br.euc <- predict(model.br.euc, newdata = br.data, type = "response")
preds.br.bcdm <- predict(model.br.bcdm, newdata = br.data, type = "response")
preds.br.bcds <- predict(model.br.bcds, newdata = br.data, type = "response")

# 3. Create and plot the ROC curve
roc_obj.br.rpd <- roc(br.data$case.control.binary, preds.br.rpd)
plot(roc_obj.br.rpd, print.auc = TRUE, col = "blue", main = "ROC for GBH",
     print.auc.x = .4,
     print.auc.y = .55,
     print.auc.cex = 1.2)

roc_obj.br.bcd <- roc(br.data$case.control.binary, preds.br.bcd)
plot(roc_obj.br.bcd, print.auc = TRUE, col = "red",
     add = T,
     print.auc.x = .4,
     print.auc.y = .47,
     print.auc.cex = 1.2)

roc_obj.br.bcdm <- roc(br.data$case.control.binary, preds.br.bcdm)
plot(roc_obj.br.bcdm, print.auc = TRUE, col = "maroon",
     add = T,
     print.auc.x = .4,
     print.auc.y = .39,
     print.auc.cex = 1.2)

roc_obj.br.bcds <- roc(br.data$case.control.binary, preds.br.bcds)
plot(roc_obj.br.bcds, print.auc = TRUE, col = "purple",
     add = T,
     print.auc.x = .4,
     print.auc.y = .31,
     print.auc.cex = 1.2)

roc_obj.br.ait <- roc(br.data$case.control.binary, preds.br.ait)
plot(roc_obj.br.ait, print.auc = TRUE, col = "yellowgreen",
     add = T,
     print.auc.x = .4,
     print.auc.y = .23,
     print.auc.cex = 1.2)

roc_obj.br.aitf <- roc(br.data$case.control.binary, preds.br.aitf)
plot(roc_obj.br.aitf, print.auc = TRUE, col = "forestgreen",
     add = T,
     print.auc.x = .4,
     print.auc.y = .15,
     print.auc.cex = 1.2)

roc_obj.br.euc <- roc(br.data$case.control.binary, preds.br.euc)
plot(roc_obj.br.euc, print.auc = TRUE, col = "orange",
     add = T,
     print.auc.x = .4,
     print.auc.y = .07,
     print.auc.cex = 1.2)

roc_obj.br.rank <- roc(br.data$case.control.binary, preds.br.rank)
plot(roc_obj.br.rank, print.auc = TRUE, col = "magenta",
     add = T,
     print.auc.x = .4,
     print.auc.y = -.01,
     print.auc.cex = 1.2)

legend("topleft", legend = c("RPD", "BCD","BCDM", "BCDS", "AIT", "AITF", "EUC", "RS"), col = c("blue", "red", "maroon", "purple", "yellowgreen", "forestgreen", "orange", "magenta"), lwd = 2, bty = "n")


####################################################################################################
####################################################################################################

####
#### Make figure showing what fraction of community is retained at different cutoffs

prop.comm.breast.mean <- vector()
prop.comm.breast.25 <- vector()
prop.comm.breast.75 <- vector()

prop.comm.agp.mean <- vector()
prop.comm.agp.25 <- vector()
prop.comm.agp.75 <- vector()

thres.vec <- seq(.4, .95, by = .01)
for(index in 1:length(thres.vec) ){
  threshold <- thres.vec[index]

  b.num.threshold <- b.num[ref, colMeans(ceiling(b.num[ref, ])) > threshold]

  prop.comm.breast.mean[index] <- mean(rowSums(b.num.threshold))
  prop.comm.breast.25[index] <- sort(rowSums(b.num.threshold))[round(.25*dim(b.num[ref, ])[1])]
  prop.comm.breast.75[index] <- sort(rowSums(b.num.threshold))[round(.75*dim(b.num[ref, ])[1])]

  abx.num.threshold <- ab.ref[, colMeans(ceiling(ab.ref)) > threshold]

  prop.comm.agp.mean[index] <- mean(rowSums(abx.num.threshold))
  prop.comm.agp.25[index] <- sort(rowSums(abx.num.threshold))[round(.25*dim(ab.ref)[1])]
  prop.comm.agp.75[index] <- sort(rowSums(abx.num.threshold))[round(.75*dim(ab.ref)[1])]
}

par(mfrow = c(1, 1))
plot(prop.comm.breast.mean ~ thres.vec,
     type = "l",
     lwd = 2,
     ylim = range(c(prop.comm.breast.25, prop.comm.breast.75, prop.comm.agp.25, prop.comm.agp.75)),
     xlab = "Reference Persistence Threshold",
     ylab = "Proportion of Community Retained",
     main = "Effect of Persistence Threshold")
lines(prop.comm.breast.25 ~ thres.vec,
      type = "l",
      lwd = 2,
      lty = 3)
lines(prop.comm.breast.75 ~ thres.vec,
      type = "l",
      lwd = 2,
      lty = 3)

## Repeat for AGP

par(mfrow = c(1, 1))
points(prop.comm.agp.mean ~ thres.vec,
       type = "l",
       lwd = 2,
       ylim = range(c(prop.comm.agp.25, prop.comm.agp.75)),
       xlab = "Reference Persistence Threshold",
       ylab = "Proportion of Community Retained",
       main = "American Gut",
       col = "red")
lines(prop.comm.agp.25 ~ thres.vec,
      type = "l",
      lwd = 2,
      lty = 3,
      col = "red")
lines(prop.comm.agp.75 ~ thres.vec,
      type = "l",
      lwd = 2,
      lty = 3,
      col = "red")

legend("bottomleft", legend = c("Breast Health", "American Gut"),
       col = c("black", "red"),
       lty = 1,
       lwd = 2,
       bty = "n")

