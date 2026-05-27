# Test/train approach for rpd vs BCD
library(pROC)

delta.auc.rpd.bcd.ab <- vector()
test.auc.rpd.ab <- vector()

delta.auc.rpd.bcd.cd <- vector()
test.auc.rpd.cd <- vector()

delta.auc.rpd.bcd.br <- vector()
test.auc.rpd.br <- vector()

train.fract <- .5

## ABX usage

set.seed(2222)

for(iter in 1:1000){

  ab.case.ind <- which(case.control.binary.abx == 1)
  ab.cont.ind <- which(case.control.binary.abx == 0)

  ab.modeltrain.ind <- c(sample(ab.case.ind, size = round(train.fract * length(ab.case.ind)) ),
                         sample(ab.cont.ind, size = round(train.fract * length(ab.cont.ind)) ) )
  ab.modeltest.ind <- c(1:length(case.control.binary.abx))[-ab.modeltrain.ind]

  # 1. Fit the GLM model
  ab.data <- as.data.frame(cbind(case.control.binary.abx, ab.test.rpd, bc.vec.ab))
  ab.data <- ab.data[ab.modeltrain.ind, ]

  model.ab.rpd <- glm(case.control.binary.abx ~ ab.test.rpd , data = ab.data, family = "binomial")
  model.ab.bcd <- glm(case.control.binary.abx ~ bc.vec.ab , data = ab.data, family = "binomial")

  summary(model.ab.rpd)
  summary(model.ab.bcd)

  coef(summary(model.ab.rpd))[, 4]
  coef(summary(model.ab.bcd))[, 4]

  ab.data.test <- ab.data[ab.modeltest.ind, ]

  # 2. Generate predicted probabilities (type = "response" is crucial)
  preds.ab.rpd <- predict(model.ab.rpd, newdata = ab.data.test, type = "response")
  preds.ab.bcd <- predict(model.ab.bcd, newdata = ab.data.test, type = "response")

  # 3. Create and plot the ROC curve
  roc_obj.ab.rpd <- roc(ab.data.test$case.control.binary.abx, preds.ab.rpd)
  roc_obj.ab.bcd <- roc(ab.data.test$case.control.binary, preds.ab.bcd)

  delta.auc.rpd.bcd.ab[iter] <- auc(roc_obj.ab.rpd) - auc(roc_obj.ab.bcd)
  test.auc.rpd.ab[iter] <- auc(roc_obj.ab.rpd)


  ## C Diff

  cd.case.ind <- which(case.control.binary.cd == 1)
  cd.cont.ind <- which(case.control.binary.cd == 0)

  cd.modeltrain.ind <- c(sample(cd.case.ind, size = round(train.fract * length(cd.case.ind))),
                         sample(cd.cont.ind, size = round(train.fract * length(cd.cont.ind))) )

  cd.modeltest.ind <- c(1:length(case.control.binary.cd))[-cd.modeltrain.ind]

  # 1. Fit the GLM model
  cd.data <- as.data.frame(cbind(case.control.binary.cd, cd.test.rpd, bc.vec.cd))
  cd.data <- cd.data[cd.modeltrain.ind, ]
  model.cd.rpd <- glm(case.control.binary.cd ~ cd.test.rpd , data = cd.data, family = "binomial")
  model.cd.bcd <- glm(case.control.binary.cd ~ bc.vec.cd , data = cd.data, family = "binomial")

  summary(model.cd.rpd)

  coef(summary(model.cd.rpd))[, 4]
  coef(summary(model.cd.bcd))[, 4]

  cd.data.test <- cd.data[cd.modeltest.ind, ]

  # 2. Generate predicted probabilities (type = "response" is crucial)
  preds.cd.rpd <- predict(model.cd.rpd, newdata = cd.data.test, type = "response")
  preds.cd.bcd <- predict(model.cd.bcd, newdata = cd.data.test, type = "response")

  # 3. Create and plot the ROC curve
  roc_obj.cd.rpd <- roc(cd.data.test$case.control.binary.cd, preds.cd.rpd)
  roc_obj.cd.bcd <- roc(cd.data.test$case.control.binary, preds.cd.bcd)

  delta.auc.rpd.bcd.cd[iter] <- auc(roc_obj.cd.rpd) - auc(roc_obj.cd.bcd)
  test.auc.rpd.cd[iter] <- auc(roc_obj.cd.rpd)

  ## Ghana Breast Health

  br.case.ind <- which(case.control.binary.gbh == 1)
  br.cont.ind <- which(case.control.binary.gbh == 0)

  br.modeltrain.ind <- c(sample(br.case.ind, size = round(train.fract * length(br.case.ind)) ),
                         sample(br.cont.ind, size = round(train.fract * length(br.cont.ind))) )
  br.modeltest.ind <- c(1:length(case.control.binary.gbh))[-br.modeltrain.ind]

  # 1. Fit the GLM model
  br.data <- as.data.frame(cbind(case.control.binary.gbh, rpd.vec.gbh, bc.vec.ref))
  br.data <- br.data[br.modeltrain.ind, ]

  model.br.rpd <- glm(case.control.binary.gbh ~ rpd.vec.gbh , data = br.data, family = "binomial")
  model.br.bcd <- glm(case.control.binary.gbh ~ bc.vec.ref , data = br.data, family = "binomial")

  summary(model.br.rpd)
  summary(model.br.bcd)

  coef(summary(model.br.rpd))[, 4]
  coef(summary(model.br.bcd))[, 4]

  br.data.test <- br.data[br.modeltest.ind, ]

  # 2. Generate predicted probabilities (type = "response" is crucial)
  preds.br.rpd <- predict(model.br.rpd, newdata = br.data.test, type = "response")
  preds.br.bcd <- predict(model.br.bcd, newdata = br.data.test, type = "response")

  # 3. Create and plot the ROC curve
  roc_obj.br.rpd <- roc(br.data.test$case.control.binary, preds.br.rpd)
  roc_obj.br.bcd <- roc(br.data.test$case.control.binary, preds.br.bcd)

  delta.auc.rpd.bcd.br[iter] <- auc(roc_obj.ab.rpd) - auc(roc_obj.br.bcd)
  test.auc.rpd.br[iter] <- auc(roc_obj.br.rpd)

}

par(mfrow = c(3, 2))

hist(test.auc.rpd.ab, breaks = seq(.5, 1, .025), main = "AGP Antibiotics", xlab = "RPD Hold-Out AUC")
mean(test.auc.rpd.ab)

hist(delta.auc.rpd.bcd.ab, breaks = seq(-.2, .25, .025), main = "AGP Antibiotics", xlab = "RPD AUC - BCD AUC")
mean(delta.auc.rpd.bcd.ab)
mean(delta.auc.rpd.bcd.ab > 0)

hist(test.auc.rpd.cd, breaks = seq(.5, 1, .025), main = "AGP C. diff", xlab = "RPD Hold-Out AUC")
mean(test.auc.rpd.cd)

hist(delta.auc.rpd.bcd.cd, breaks = seq(-.2, .25, .025), main = "AGP C. diff", xlab = "RPD AUC - BCD AUC")
mean(delta.auc.rpd.bcd.cd)
mean(delta.auc.rpd.bcd.cd > 0)

hist(test.auc.rpd.br, breaks = seq(.5, 1, .025), main = "GBH", xlab = "RPD Hold-Out AUC")
mean(test.auc.rpd.br)

hist(delta.auc.rpd.bcd.br, breaks = seq(-.2, .25, .025), main = "GBH", xlab = "RPD AUC - BCD AUC")
mean(delta.auc.rpd.bcd.br)
mean(delta.auc.rpd.bcd.br > 0)


