import numpy as np

def ref_rank(val, vec):
    """
    Calculate where a value falls within another vector (proportion of values in vec less than val).
    
    Parameters
    ----------
    val : float
        Value to rank
    vec : array-like
        Reference vector
        
    Returns
    -------
    float
        Proportion of values in vec less than val
    """
    vec_finite = vec[np.isfinite(vec)]
    if len(vec_finite) == 0:
        return np.nan
    return np.sum(val > vec_finite) / len(vec_finite)

def ref_rank_self(val, vec):
    """
    Calculate where a value falls within another vector, excluding one matching value.
    
    Parameters
    ----------
    val : float
        Value to rank
    vec : array-like
        Reference vector
        
    Returns
    -------
    float
        Proportion of values in vec (excluding one instance of val) less than val
    """
    # Find which value in vec is equal to the value passed. If there are multiple, just take 1
    remove_idx = np.where(vec == val)[0]
    if len(remove_idx) > 0:
        # Remove the first matching value from the reference vector
        vec = np.delete(vec, remove_idx[0])
    # Proceed as before
    vec_finite = vec[np.isfinite(vec)]
    if len(vec_finite) == 0:
        return np.nan
    return np.sum(val > vec_finite) / len(vec_finite)

def calc_rpd(ref, test, ref_pers, taxon_scores=False):
    """
    Calculate Ratio Percentile Deviation (RPD) for each test sample, given a reference dataset.
    
    Parameters
    ----------
    ref : array-like, shape (n_ref_samples, n_taxa)
        Reference dataset (rows=samples, columns=taxa)
    test : array-like, shape (n_test_samples, n_taxa)
        Test dataset (rows=samples, columns=taxa)
    ref_pers : float
        Reference persistence threshold (proportion of samples where taxon must be present)
    taxon_scores : bool, optional (default=False)
        If True, return both taxon-level and sample-level RPD scores
        
    Returns
    -------
    If taxon_scores is False:
        sample_rpd : array, shape (n_test_samples,)
            RPD score for each test sample
    If taxon_scores is True:
        tuple (taxon_rpd, sample_rpd)
            taxon_rpd : array, shape (n_test_samples, n_taxa)
                Taxon-level RPD deviations
            sample_rpd : array, shape (n_test_samples,)
                Sample-level RPD scores (mean of taxon-level deviations)
                
    Notes
    -----
    This function implements the RPD method as described in the original R implementation.
    The algorithm:
    1. Filters taxa by persistence in the reference dataset
    2. For each taxon, calculates ratios to all other taxa in both reference and test sets
    3. Determines percentile ranks of test ratios within reference ratio distributions
    4. Computes average deviation from median (0.5) as the RPD score
    """
    # Convert inputs to numpy arrays
    ref = np.asarray(ref)
    test = np.asarray(test)
    
    # Check that ref and test have same number of taxa
    if ref.shape[1] != test.shape[1]:
        raise ValueError("ref and test datasets must contain the same taxa")
    
    # Check whether ref and test are the same, for use in choosing ref.rank function
    self_ref = np.array_equal(ref, test)
    
    # Subset reference dataset to only taxa with the requisite persistence
    # Persistence = proportion of samples where taxon is present (abundance > 0)
    persistence = np.mean(ref > 0, axis=0)
    keep_taxa = persistence > ref_pers
    ref_sub = ref[:, keep_taxa]
    
    # Match up the test dataset to the taxa retained in the ref.sub dataset
    test_sub = test[:, keep_taxa]
    
    # Initialize output matrix for taxon-level rank variances
    n_test_samples, n_kept_taxa = test_sub.shape
    ratio_diff_mat = np.full((n_test_samples, n_kept_taxa), np.nan)
    
    # Process each taxon (column)
    for col_idx in range(n_kept_taxa):
        # Pull out abundance of focal taxon in ref dataset
        focal_col_ref = ref_sub[:, col_idx]
        
        # Create ratios of other taxa to focal taxon within ref dataset
        # Avoid division by zero: where focal_col_ref is 0, set ratio to inf
        with np.errstate(divide='ignore', invalid='ignore'):
            tax_rat_ref = ref_sub / focal_col_ref[:, np.newaxis]
            # Replace inf/nan from division by zero with nan for ranking
            tax_rat_ref[~np.isfinite(tax_rat_ref)] = np.nan
        
        # Pull out abundance of taxon in test dataset
        focal_tax_test = test_sub[:, col_idx]
        
        # Calculate ratios of focal taxon within test dataset
        with np.errstate(divide='ignore', invalid='ignore'):
            tax_rat_test = test_sub / focal_tax_test[:, np.newaxis]
            tax_rat_test[~np.isfinite(tax_rat_test)] = np.nan
        
        # Create a matrix to hold ranks for this focal taxon
        tax_rat_rank_mat = np.full((n_test_samples, n_kept_taxa), np.nan)
        
        if not self_ref:
            # Write a loop to fill in taxon level percentile using ref.rank
            for taxon_idx in range(n_kept_taxa):
                # Extract vectors for this taxon pair
                test_vec = tax_rat_test[:, taxon_idx]
                ref_vec = tax_rat_ref[:, taxon_idx]
                
                # Apply ref.rank to each value in test_vec
                ranks = np.array([ref_rank(val, ref_vec) for val in test_vec])
                tax_rat_rank_mat[:, taxon_idx] = ranks
            
            # Remove values for focal taxon (set to NaN)
            tax_rat_rank_mat[:, col_idx] = np.nan
        else:
            # Same loop but with ref.rank.self
            for taxon_idx in range(n_kept_taxa):
                test_vec = tax_rat_test[:, taxon_idx]
                ref_vec = tax_rat_ref[:, taxon_idx]
                
                ranks = np.array([ref_rank_self(val, ref_vec) for val in test_vec])
                tax_rat_rank_mat[:, taxon_idx] = ranks
            
            # Remove values for focal taxon
            tax_rat_rank_mat[:, col_idx] = np.nan
        
        # Calculate average difference from median percentile for this taxon
        # Only consider test samples where the focal taxon is present (>0)
        present_mask = test_sub[:, col_idx] > 0
        if np.any(present_mask):
            # For each sample, compute mean absolute deviation from 0.5 across taxa
            dev_from_median = np.abs(tax_rat_rank_mat[present_mask, :] - 0.5)
            # Mean across taxa (ignoring NaN values)
            mean_dev = np.nanmean(dev_from_median, axis=1)
            ratio_diff_mat[present_mask, col_idx] = mean_dev
    
    # Calculate the mean of the average deviations in each sample
    rpd_vector = np.nanmean(ratio_diff_mat, axis=1)
    
    # Replace any values from samples with no overlapping taxa with 0.5
    zero_taxa_mask = np.sum(test_sub > 0, axis=1) == 0
    rpd_vector[zero_taxa_mask] = 0.5
    
    if taxon_scores:
        return ratio_diff_mat, rpd_vector
    return rpd_vector