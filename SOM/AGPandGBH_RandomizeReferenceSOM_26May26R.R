# Repeat main text analyses 100 times with randomized references

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

#### AGP setup

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



####
#### Results loop

runs <- 100

br.rpd.auc <- vector()
br.bcd.auc <- vector()
ab.rpd.auc <- vector()
ab.bcd.auc <- vector()
cd.rpd.auc <- vector()
cd.bcd.auc <- vector()

set.seed(3333)

for(run in 1:runs){

####
#### GBH

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

# # Calculate rpd for controls
# control.rpd.gbh <- calc.rpd(ref = b.num[ref, ], test = b.num[control, ], ref.pers = ref.threshold)
#
# # Calculate rpd for cases
# case.rpd.gbh <- calc.rpd(ref = b.num[ref, ], test = b.num[case, ], ref.pers = ref.threshold)

rpd.vec.gbh <- calc.rpd(ref = b.num[ref, ], test = b.num[c(case, control), ], ref.pers = ref.threshold)

case.control.vec.gbh <- c(rep("case", times = length(case)), rep("control", times = length(control)))

table(case.control.vec.gbh)

# # Calculate summary stats on rpd values
# mean(case.rpd.gbh)
# sd(case.rpd.gbh)
#
# mean(control.rpd.gbh)
# sd(control.rpd.gbh)

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

# # Implement same logistic regression
# summary(glm(as.factor(case.control.binary.gbh) ~ bc.vec.ref, family = "binomial"))
#
# ## Calculate diversity metrics
# ref.div.gbh <- vegan::diversity(b.num[ref, ])
# control.div.gbh <- vegan::diversity(b.num[control, ])
# case.div.gbh <- vegan::diversity(b.num[case, ])
# case.control.div.gbh <- c(case.div.gbh, control.div.gbh)

####
#### AGP data

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
# summary(glm(as.factor(case.control.binary.abx) ~ ab.test.rpd, family = "binomial"))

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
cd.ref.mean <- apply(cd.ref, 2, mean)

bc.vec.cd <- as.matrix(vegdist(rbind(cd.ref.mean, cd.test)))[-1, 1]

# Implement logistic regression
# summary(glm(as.factor(case.control.binary.cd) ~ bc.vec.cd, family = "binomial"))

####################################################################################################
####################################################################################################

# Make figures of rpd versus diversity and BCD

# First, do a check to make sure that sequencing depth does not imapct diversity too much:
# Calculate raw diversity, then remove all values below the threshold of the smallest value consistently observed in the data, and re-calculate
# If correlation of diversity metrics in the two scenarios is high, then rare taxa do not affect diversity

# d.div <- vegan::diversity(d.sub1)
#
# # Plot min detection in each sample
# min.non.zero <- function(vec){
#   min.pos <- min(vec[vec>0])
# }
#
# # Find detection threshold by looking at the smallest observed values in each sample
# hist(log(apply(d, 1, min.non.zero)))
# hist(log(apply(d[cd.test.indices,], 1, min.non.zero)))
# exp(-8) # this value is pretty consistently detected
#
# # Create a dataset removing all values smaller than e^-8
# d.cut8 <- d.sub1
# d.cut8[d.cut8 < exp(-8)] <- 0
#
# # Calculate diversity on this thresholded dataset
# d8.div <- vegan::diversity(d.cut8)
#
# # Look at correlation of diversity calculated on two datasets
# plot(d.div, d8.div) # excellent, cutting out small values has almost no impact
# cor(d.div, d8.div) # 0.998. Valid to use raw diversity
#
# ####
# #### Make RPD vs diversity and RPD vs BCD figures:
#
# par(mfrow = c(3, 2))
# plot(ab.test.rpd ~ d.div[ab.test.indices],
#      ylim = range(c(ab.test.rpd, cd.test.rpd, rpd.vec.gbh + .02), na.rm =  T),
#      xlim = range(c(d.div, case.control.div.gbh)),
#      col = c(rep(adjustcolor("royalblue1", alpha = .5), times = sum(case.control.binary.abx == 0) ), rep(adjustcolor("royalblue4", alpha = .5), times = sum(case.control.binary.abx == 1) ) ),
#      xlab = "Diversity",
#      ylab = "Sample RPD",
#      main = "AGP Antibiotics",
#      pch = 16,
#      cex= 1.2)
# text(labels = paste("r = ", round(cor(ab.test.rpd , d.div[ab.test.indices], use = "pairwise.complete.obs"), digits = 3) ) , x = 4.5, y = .48, cex = 1.2)
#
# legend("bottomleft", legend = c("case", "control"),
#        col = c(adjustcolor("royalblue4", alpha = .5), adjustcolor("royalblue1", alpha = .5)),
#        pch = 16,
#        cex = 1.2)
#
# cor(ab.test.rpd , d.div[ab.test.indices], use = "pairwise.complete.obs")
#
# plot(ab.test.rpd ~ bc.vec.ab,
#      xlim = range(c(bc.vec.ab, bc.vec.cd, bc.vec.ref), na.rm =  T),
#      ylim = range(c(ab.test.rpd, cd.test.rpd, rpd.vec.gbh + .02), na.rm =  T),
#      col = c(rep(adjustcolor("royalblue1", alpha = .5), times = sum(case.control.binary.abx == 0) ), rep(adjustcolor("royalblue4", alpha = .5), times = sum(case.control.binary.abx == 1) ) ),
#      xlab = "BCD to Reference Centroid",
#      ylab = "Sample RPD",
#      main = "AGP Antibiotics",
#      pch = 16,
#      cex= 1.2)
# text(labels = paste("r = ", formatC(signif(cor(ab.test.rpd , bc.vec.ab, use = "pairwise.complete.obs"), digits = 3), digits = 3, format="fg", flag="#" ) ) , x = .47, y = .48, cex = 1.2)
#
# cor(ab.test.rpd , bc.vec.ab, use = "pairwise.complete.obs")
#
# plot(cd.test.rpd ~ d.div[cd.test.indices],
#      ylim = range(c(ab.test.rpd, cd.test.rpd, rpd.vec.gbh + .02), na.rm =  T),
#      xlim = range(c(d.div, case.control.div.gbh)),
#      col = c(rep(adjustcolor("red1", alpha = .4), times = sum(case.control.binary.cd == 0) ), rep(adjustcolor("red4", alpha = .4), times = sum(case.control.binary.cd == 1) ) ),
#      xlab = "Diversity",
#      ylab = "Sample RPD",
#      main = "AGP C. diff",
#      pch = 16,
#      cex= 1.2)
# text(labels = paste("r = ", round(cor(cd.test.rpd , d.div[cd.test.indices], use = "pairwise.complete.obs"), digits = 3) ) , x = 4.5, y = .48, cex = 1.2)
# legend("bottomleft", legend = c("case", "control"),
#        col = c(adjustcolor("red4", alpha = .4), adjustcolor("red1", alpha = .4)),
#        pch = 16,
#        cex = 1.2)
#
# cor(cd.test.rpd , d.div[cd.test.indices], use = "pairwise.complete.obs")
#
# plot(cd.test.rpd ~ bc.vec.cd,
#      xlim = range(c(bc.vec.ab, bc.vec.cd, bc.vec.ref), na.rm =  T),
#      ylim = range(c(ab.test.rpd, cd.test.rpd, rpd.vec.gbh + .02), na.rm =  T),
#      col = c(rep(adjustcolor("red1", alpha = .4), times = sum(case.control.binary.cd == 0) ), rep(adjustcolor("red4", alpha = .4), times = sum(case.control.binary.cd == 1) ) ),
#      xlab = "BCD to Reference Centroid",
#      ylab = "Sample RPD",
#      main = "AGP C. diff",
#      pch = 16,
#      cex= 1.2)
# text(labels = paste("r = ", formatC(round(cor(cd.test.rpd , bc.vec.cd, use = "pairwise.complete.obs"), digits = 3), digits=3, format="fg", flag="#")) , x = .47, y = .48, cex = 1.2)
#
# cor(cd.test.rpd , bc.vec.cd, use = "pairwise.complete.obs")
#
# plot(rpd.vec.gbh ~ case.control.div.gbh,
#      ylim = range(c(ab.test.rpd, cd.test.rpd, rpd.vec.gbh + .02), na.rm =  T),
#      xlim = range(c(d.div, case.control.div.gbh)),
#      col = c(rep(adjustcolor("darkgoldenrod4", alpha = .4), times = sum(case.control.vec.gbh == "case") ), rep(adjustcolor("darkgoldenrod1", alpha = .4), times = sum(case.control.vec.gbh == "control") ) ),
#      xlab = "Diversity",
#      ylab = "Sample RPD",
#      main = "Ghana Breast Health",
#      pch = 16,
#      cex= 1.2)
# text(labels = paste("r = ", round(cor(rpd.vec.gbh , case.control.div.gbh, use = "pairwise.complete.obs"), digits = 3) ) , x = 4.5, y = .48, cex = 1.2)
#
# legend("bottomleft", legend = c("case", "control"),
#        col = c(adjustcolor("darkgoldenrod4", alpha = .4), adjustcolor("darkgoldenrod1", alpha = .4)),
#        pch = 16,
#        cex = 1.2)
#
# cor(rpd.vec.gbh , case.control.div.gbh, use = "pairwise.complete.obs")
#
# plot(rpd.vec.gbh ~ bc.vec.ref,
#      xlim = range(c(bc.vec.ab, bc.vec.cd, bc.vec.ref), na.rm =  T),
#      ylim = range(c(ab.test.rpd, cd.test.rpd, rpd.vec.gbh + .02), na.rm =  T),
#      col = c(rep(adjustcolor("darkgoldenrod4", alpha = .4), times = sum(case.control.vec.gbh == "case") ), rep(adjustcolor("darkgoldenrod1", alpha = .4), times = sum(case.control.vec.gbh == "control") ) ),
#      xlab = "BCD to Reference Centroid",
#      ylab = "Sample RPD",
#      main = "Ghana Breast Health",
#      pch = 16,
#      cex= 1.2)
# text(labels = paste("r = ", round(cor(rpd.vec.gbh , bc.vec.ref, use = "pairwise.complete.obs"), digits = 3) ) , x = .47, y = .48, cex = 1.2)
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
# par(cex.lab = 1.4, cex.main = 1.5)
#
# vioplot(ab.test.rpd ~ case.control.binary.abx,
#         names = c("None", "Recent Antibiotics"),
#         xlab = "Antibiotic Use",
#         ylab = "Sample RPD",
#         col = c(adjustcolor("royalblue1", alpha = .5), adjustcolor("royalblue4", alpha = .5)) ,
#         cex.axis = 1.3,
#         cex.lab = 1.4,
#         cex.main = 1.6,
#         main = "AGP Antibiotics",
#         ylim = range(c(ab.test.rpd, cd.test.rpd, rpd.vec.gbh)))
#
# vioplot(cd.test.rpd ~ meta.cd.test,
#         names = c( "No diagnosis", "Diagnosed with C. diff"),
#         xlab = "C. diff Diagnosis Status",
#         ylab = "Sample RPD",
#         cex.axis = 1.3,
#         cex.lab = 1.4,
#         cex.main = 1.6,
#         main = "AGP C. diff",
#         col = c(adjustcolor("red1", alpha = .5), adjustcolor("red4", alpha = .5)) ,
#         ylim = range(c(ab.test.rpd, cd.test.rpd, rpd.vec.gbh)))
#
# vioplot(rpd.vec.gbh ~ case.control.vec.gbh,
#         names = c("Healthy", "Breast Disease"),
#         xlab = "Breast Disease Status",
#         ylab = "Sample RPD",
#         main = "Ghana Breast Health",
#         cex.axis = 1.3,
#         cex.lab = 1.4,
#         cex.main = 1.6,
#         col = c(adjustcolor("darkgoldenrod1", alpha = .5), adjustcolor("darkgoldenrod4", alpha = .5)) ,
#         ylim = range(c(ab.test.rpd, cd.test.rpd, rpd.vec.gbh)))
#

#### ROC Curves

## ABX usage

# 1. Fit the GLM model
ab.data <- as.data.frame(cbind(case.control.binary.abx, ab.test.rpd, bc.vec.ab))
model.ab.rpd <- glm(case.control.binary.abx ~ ab.test.rpd , data = ab.data, family = "binomial")
model.ab.bcd <- glm(case.control.binary.abx ~ bc.vec.ab , data = ab.data, family = "binomial")

# 2. Generate predicted probabilities (type = "response" is crucial)
preds.ab.rpd <- predict(model.ab.rpd, newdata = ab.data, type = "response")
preds.ab.bcd <- predict(model.ab.bcd, newdata = ab.data, type = "response")

# 3. Create and plot the ROC curve
roc_obj.ab.rpd <- roc(ab.data$case.control.binary.abx, preds.ab.rpd)

roc_obj.ab.bcd <- roc(ab.data$case.control.binary, preds.ab.bcd)


## C Diff

# 1. Fit the GLM model
cd.data <- as.data.frame(cbind(case.control.binary.cd, cd.test.rpd, bc.vec.cd))
model.cd.rpd <- glm(case.control.binary.cd ~ cd.test.rpd , data = cd.data, family = "binomial")
model.cd.bcd <- glm(case.control.binary.cd ~ bc.vec.cd , data = cd.data, family = "binomial")

# 2. Generate predicted probabilities (type = "response" is crucial)
preds.cd.rpd <- predict(model.cd.rpd, newdata = cd.data, type = "response")
preds.cd.bcd <- predict(model.cd.bcd, newdata = cd.data, type = "response")

# 3. Create and plot the ROC curve
roc_obj.cd.rpd <- roc(cd.data$case.control.binary.cd, preds.cd.rpd)

roc_obj.cd.bcd <- roc(cd.data$case.control.binary, preds.cd.bcd)


## Breast Health

# 1. Fit the GLM model
br.data <- as.data.frame(cbind(case.control.binary.gbh, rpd.vec.gbh, bc.vec.ref))
model.br.rpd <- glm(case.control.binary.gbh ~ rpd.vec.gbh , data = br.data, family = "binomial")
model.br.bcd <- glm(case.control.binary.gbh ~ bc.vec.ref , data = br.data, family = "binomial")

# 2. Generate predicted probabilities (type = "response" is crucial)
preds.br.rpd <- predict(model.br.rpd, newdata = br.data, type = "response")
preds.br.bcd <- predict(model.br.bcd, newdata = br.data, type = "response")

# 3. Create and plot the ROC curve
roc_obj.br.rpd <- roc(br.data$case.control.binary, preds.br.rpd)
roc_obj.br.bcd <- roc(br.data$case.control.binary, preds.br.bcd)


ab.rpd.auc[run] <- auc(roc_obj.ab.rpd)
ab.bcd.auc[run] <- auc(roc_obj.ab.bcd)
cd.rpd.auc[run] <- auc(roc_obj.cd.rpd)
cd.bcd.auc[run] <- auc(roc_obj.cd.bcd)
br.rpd.auc[run] <- auc(roc_obj.br.rpd)
br.bcd.auc[run] <- auc(roc_obj.br.bcd)

if(run %% 5 == 0) print(run)
}

par(mfrow = c(1, 3))
hist(ab.rpd.auc - ab.bcd.auc, xlab = "RPD AUC - BCD AUC", main = "AGP Antibiotics", breaks = 6)

hist(cd.rpd.auc - cd.bcd.auc, xlab = "RPD AUC - BCD AUC", main = "AGP C. diff", breaks = 6)

hist(br.rpd.auc - br.bcd.auc, xlab = "RPD AUC - BCD AUC", main = "GBH", breaks = 10)

mean(ab.rpd.auc - ab.bcd.auc)
mean(cd.rpd.auc - cd.bcd.auc)
mean(br.rpd.auc - br.bcd.auc)

min(ab.rpd.auc - ab.bcd.auc)
min(cd.rpd.auc - cd.bcd.auc)
min(br.rpd.auc - br.bcd.auc)

