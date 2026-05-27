## Lake Mendota RPD Analysis

# Load packages and data
library(lubridate)
library(vegan)
library(limony)
data("limony")

# Save OTU table, taxonomy, and sample dates
otu.table <- limony$av$seqID / 100
otu.taxonomy <- limony$names$seqID
sample.dates <- convert.sample.names.to.dates(sample.names = colnames(otu.table))

# Reformat OTU table
tab <- t(otu.table)

####
#### Read in dissolved nutrient data and do quality control
np <- read.csv("~/Desktop/MethodsRRV/MEnutrients.csv", header = T)
# Subset just to consistently sampled months
start.month <- 5
end.month <- 10
np <- np[np$Month >= start.month & np$Month <= end.month, ]
# Remove any replicates not A or B
np <- np[np$Biol.Rep %in% c("A", "B"), ]
# Remove any non standard depths
np <- np[!grepl("D", np$Sample.Name), ]
# Create a date entry in np
np$date <- paste(paste(np$Year, np$Month, sep = "-"), np$Day, sep = "-")
np$day.of.year <- yday(as.Date(np$date, format = "%Y-%m-%d"))

## Identify outliers by looking at SD of log transformed replicate nutrient samples

par(mfrow = c(1, 1))
hist(tapply(log(np$TDN.ug.L), np$Sample.Name, sd) , breaks = 25)
hist(tapply(log(np$TDP.ug.L), np$Sample.Name, sd) , breaks = 25)

outlier.n <- which(tapply(log(np$TDN.ug.L), np$Sample.Name, sd) > .25)
outlier.p <- which(tapply(log(np$TDP.ug.L), np$Sample.Name, sd) > .6)

# Inspect these data
np[np$Sample.Name %in% c(names(outlier.n), names(outlier.p)), ]

## Remove outlier points
# Take out p concentrations lower than 3
sum(np$TDP.ug.L < 3)
np$TDP.ug.L[np$TDP.ug.L < 3 ] <- NA

# Take out n concentrations lower than 400
sum(np$TDN.ug.L < 400)
np$TDN.ug.L[np$TDN.ug.L < 400 ] <- NA

# Manually remove remaining anomalous data points
np$TDN.ug.L[np$Sample.Name ==  "ME23Aug2013" & np$Biol.Rep == "B" ] <- NA
np$TDP.ug.L[np$Sample.Name ==  "ME23Aug2013" & np$Biol.Rep == "B" ] <- NA
np$TDN.ug.L[np$Sample.Name ==  "ME22Sep2017" & np$Biol.Rep == "A" ] <- NA
np$TDN.ug.L[np$Sample.Name ==  "ME06Jul2016" & np$Biol.Rep == "A" ] <- NA

####
#### Format community composition data

    # Set parameters for start and end years for analysis
    start.year <- 2013
    end.year <- 2018

    # Set ref persistence value for RPD
    ref.pers <- .85 # .85 and .9 both good

    # Create subsetted dataset using year and month cutoffs
    tab.sub <- tab[year(sample.dates) >= start.year & year(sample.dates) <= end.year & month(sample.dates) >= start.month & month(sample.dates) <= end.month, ]

    dim(tab.sub)

    # Calculate RPD based on OTU table
    otu.rpd <- calc.rpd(ref = tab.sub, test = tab.sub, ref.pers = ref.pers)

    # Find how many taxa are retained
    sum(colMeans(ceiling(tab.sub)) > ref.pers)

    # Save the dates associated with the samples
    rpd.dates <- as.Date(convert.sample.names.to.dates(sample.names = rownames(tab.sub)))

####
#### Create metrics for divergence of N and P data

    # Subset np to only years analyzed. Already was subset to months analyzed.

    ref.n.mean <- mean(log(np$TDN.ug.L[np$Year %in% start.year:end.year ]), na.rm = T)
    ref.n.sd <- sd(log(np$TDN.ug.L[np$Year %in% start.year:end.year ]), na.rm = T)

    ref.p.mean <- mean(log(np$TDP.ug.L[np$Year %in% start.year:end.year ]), na.rm = T)
    ref.p.sd <- sd(log(np$TDP.ug.L[np$Year %in% start.year:end.year ]), na.rm = T)

    ## Create z scores for test data

    # subset data to test year
    test.np <- np[np$Year %in% start.year:end.year, ]

    # Average n and p data per date
    test.n.agg <- aggregate(TDN.ug.L ~ Sample.Name, data = test.np, FUN = mean, na.rm = T)
    test.p.agg <- aggregate(TDP.ug.L ~ Sample.Name, data = test.np, FUN = mean, na.rm = T)

    # Re-order dates to be in chronological order
    test.n.ord <- test.n.agg[match(unique(test.np$Sample.Name), test.n.agg$Sample.Name), ]
    test.p.ord <- test.p.agg[match(unique(test.np$Sample.Name), test.p.agg$Sample.Name), ]

    # Create z score for each log-transformed date
    test.z.n <- ( log(test.n.ord$TDN.ug.L) - ref.n.mean ) / ref.n.sd
    test.z.p <- ( log(test.p.ord$TDP.ug.L) - ref.p.mean ) / ref.p.sd

    # Create date vector for each
    n.dates.test <- as.Date(convert.sample.names.to.dates(sample.names = substring(test.n.ord$Sample.Name, first = 3, last = 11)))
    p.dates.test <- as.Date(convert.sample.names.to.dates(sample.names = substring(test.p.ord$Sample.Name, first = 3, last = 11)))

####
#### Match up N and P values with RPD values by pairing same dates

      ## N
      # Find dates of N in rpd
      match.rpd.n <- match(rpd.dates, n.dates.test)
      n.in.rpd <- test.z.n[match.rpd.n]

      final.n.dates <- n.dates.test[match.rpd.n]

      ## P
      # Find dates of P in rpd
      match.rpd.p <- match(rpd.dates, p.dates.test)
      p.in.rpd <- test.z.p[match.rpd.p]

      final.p.dates <- p.dates.test[match.rpd.p]

####
#### Calculate BCD to conduct comparable analysis

      ## Calculate distance from centroid to run parallel analysis

      # Take the centroid of the full dataset
      mean.tab.p.bc <- apply(tab.sub, 2, mean)

      # Find communities with matching nutrient samples
      tab.p.bc <- tab.sub[ as.Date(convert.sample.names.to.dates(rownames(tab.sub))) %in% final.p.dates, ]

      # Take the vector of BCD from each sample to the centroid
      bc.tab.p <- as.matrix(vegdist(rbind(mean.tab.p.bc, tab.p.bc)))[1, -1]


      ## Compare variance explained in the multiple regressions
      summary(lm(otu.rpd ~ abs(n.in.rpd) + abs(p.in.rpd)))
      summary(lm(bc.tab.p ~ na.omit(abs(n.in.rpd)) + na.omit(abs(p.in.rpd))))


####
#### Make figure illustrating multiple regression

      num.p.cols <- 200
      scaled.p.rpd.col <- ceiling(na.omit(abs(p.in.rpd)) * num.p.cols / max(na.omit(abs(p.in.rpd))))

      ## Create figures to show results
      par(mfrow = c(1, 2), cex.axis = 1.15, cex.lab = 1.25, cex.main = 1.25)
      plot(abs(n.in.rpd)[!is.na(n.in.rpd)] , otu.rpd[!is.na(n.in.rpd)],
           xlab = "Dissolved N Deviation",
           ylab = "Sample RPD",
           main = "RPD vs. Nutrient Variability",
           #col = ((hcl.colors(200, "YlOrRd")))[scaled.p.rpd.col ],
           col = adjustcolor(heat.colors(210)[scaled.p.rpd.col ], alpha.f = 1),
           pch = 16)
      abline(lm(otu.rpd[!is.na(n.in.rpd)] ~ abs(n.in.rpd)[!is.na(n.in.rpd)]),
             lty = 3, lwd = 3)
      text(1.4, .1885, labels = expression(R^2 == 0.33), cex = 1.3)

      plot(abs(na.omit(n.in.rpd)) , bc.tab.p,
           xlab = "Dissolved N Deviation",
           ylab = "BCD to Dataset Centroid",
           main = "BCD vs. Nutrient Variability",
           #col = ((hcl.colors(200, "YlOrRd")))[scaled.p.rpd.col ],
           col = adjustcolor(heat.colors(210)[scaled.p.rpd.col ], alpha.f = 1),
           pch = 16)
      abline(lm(bc.tab.p ~ abs(na.omit(n.in.rpd)) ),
             lty = 3, lwd = 3)
      text(1.4, .398, labels = expression(R^2 == 0.17), cex = 1.3)

      cor(otu.rpd[!is.na(p.in.rpd)], bc.tab.p)

####
#### Make descriptive figure for N and P and microbial community composition

## Plot of seasonal nitrogen and phosphorus dynamics with log-transformed scale

# Use date entry to make seasonal plot with different colors for each year
par(mfrow = c(3,1), mar = c(4, 4, 2, 2))
plot(np$TDN.ug.L ~ np$day.of.year,
     col = as.factor(np$Year),
     pch = 16,
     log = "y",
     xlab = "Day of Year",
     ylab = "Total Dissolved N",
     cex.axis = 1.3, cex.lab = 1.3)
legend("bottomleft", legend = unique(np$Year),
       col = as.factor(unique(np$Year)),
       pch = 16,
       bty = "n")

plot(np$TDP.ug.L ~ np$day.of.year,
     col = as.factor(np$Year),
     pch = 16,
     log = "y",
     xlab = "Day of Year",
     ylab = "Total Dissolved P",
     cex.axis = 1.3, cex.lab = 1.3)


## Plot PCA of samples

par(mfrow = c(1,1))
dim(tab.sub)
tab.pca <- tab.sub[ , colMeans(ceiling(tab.sub) ) > .01]
dim(tab.pca)

me.pca <- prcomp((tab.pca))

pca.scores <- me.pca$x[,1:2]
plot(pca.scores[,1], pca.scores[,2],
    col = "white",
    xlab = "PC1, 32.6% variability",
    ylab = "PC2, 13.6% variability",
    cex.lab = 1, cex.axis = 1)
text( pca.scores[,1], pca.scores[,2],
      col = as.numeric(as.factor(substr(rownames(tab.pca), 6, 9))),
      labels = substr(rownames(tab.pca), 3, 5),
      lwd = 1.5,
      font = 2)

summary(me.pca)$importance[, 1:2]


