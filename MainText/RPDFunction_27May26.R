#### Create a function to calculate RPD for each test sample, given a reference dataset

## Create supporting functions

# Create function to use with lapply that calculates where a position falls within another vector
ref.rank <- function(val, vec){
  vec.finite <- vec[is.finite(vec)] # Remove NaN and Inf values
  ref.rank.fract <- sum(val > vec.finite) / length(vec.finite)
  return(ref.rank.fract)
}

# Modify above function to use when ref and test datasets are the same
ref.rank.self <- function(val, vec){
  # Find which value in vec is equal to the value passed. If there are multiple, just take 1
  remove.val <- which(vec == val)[1]
  # Remove the value from the reference vector
  vec <- vec[-remove.val]
  # Proceed as before
  vec.finite <- vec[is.finite(vec)] # Remove NaN and Inf values
  ref.rank.fract <- sum(val > vec.finite) / length(vec.finite)
  return(ref.rank.fract)
}

## Create main RPD function

# Change optional argument taxon.scores to be T if wanting taxon-level RPD as well

calc.rpd <- function(ref, test, ref.pers, taxon.scores = F){

  # Check that ref and test have same number of taxa
  if (dim(ref)[2] != dim(test)[2]) {
    stop("ref and test datasets must contain the same taxa")
  }

  # Check whether ref and test are the same, for use in choosing ref.rank function
  if(identical(ref, test)) {
    self.ref <- T
  } else {
      self.ref <- F
      }

  # Subset reference dataset to only taxa with the requisite persistence
  ref.sub <- ref[, colMeans(ref > 0) > ref.pers]

  # Match up the test dataset to the taxa retained in the ref.sub dataset
  test.sub <- test[, colMeans(ref > 0) > ref.pers]

  # Create taxon-level rank variances.
  # The loop will store the average difference from median of ratio percentiles for each taxon for each observation in the test dataset
  ratio.diff.mat <- matrix(nrow = dim(test.sub)[1], ncol = dim(test.sub)[2])
  colnames(ratio.diff.mat) <- colnames(test.sub)

  for(col in 1:dim(ref.sub)[2]){
    # Pull out abundance of focal taxon in ref dataset
    focal.col.ref <- ref.sub[, col]

    # Create ratios of other taxa to focal taxon within ref dataset
    tax.rat.ref <- sweep(ref.sub, 1, focal.col.ref, FUN = '/')

    # Pull out abundance of taxon in test dataset
    focal.tax.test <- test.sub[, col]

    # Calculate ratios of focal taxon within test dataset
    tax.rat.test <- sweep(test.sub, 1, focal.tax.test, FUN = '/')

    # Create a matrix to hold ranks
    tax.rat.rank.mat <- matrix(nrow = dim(test.sub)[1], ncol = dim(test.sub)[2])

    if(self.ref == F) {
      # Write a for loop to fill in taxon level percentile
      for(taxon in 1:dim(tax.rat.rank.mat)[2]){
        tax.rat.rank.mat[, taxon] <- unlist(lapply(tax.rat.test[, taxon], ref.rank, vec = tax.rat.ref[, taxon]))
      }
      # Remove values for focal taxon, which will always be the same
      tax.rat.rank.mat[, col] <- NA
      } else {
        #Same loop but with ref.rank.self
        for(taxon in 1:dim(tax.rat.rank.mat)[2]){
          tax.rat.rank.mat[, taxon] <- unlist(lapply(tax.rat.test[, taxon], ref.rank.self, vec = tax.rat.ref[, taxon]))
        }
        # Remove values for focal taxon, which will always be the same
        tax.rat.rank.mat[, col] <- NA

      }

  # Calculate average difference from median percentile
  ratio.diff.mat[test.sub[, col] > 0, col] <- rowMeans(abs(tax.rat.rank.mat[test.sub[, col] > 0, ] - .5), na.rm = T)

  }

  # Calculate the mean of the average deviations in each sample
  rpd.vector <- apply(ratio.diff.mat, 1, mean, na.rm = T)

  # Replace any values from samples with no overlapping taxa with 0.5
  rpd.vector[which(rowSums(test.sub > 0) == 0)] <- 0.5

  if(taxon.scores == T) {
    results.list <- list(taxon.rpd = ratio.diff.mat, sample.rpd = rpd.vector)
    return(results.list) } else {
      return(rpd.vector)
    }

}

