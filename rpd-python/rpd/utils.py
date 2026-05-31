"""
Utility functions for RPD package
"""

import numpy as np
import pandas as pd

def validate_input_data(ref, test):
    """
    Validate input data for RPD calculation.
    
    Parameters
    ----------
    ref : array-like
        Reference dataset
    test : array-like
        Test dataset
        
    Returns
    -------
    tuple
        (ref_array, test_array) as numpy arrays
        
    Raises
    ------
    ValueError
        If inputs are invalid
    """
    # Convert to numpy arrays if needed
    ref_array = np.asarray(ref)
    test_array = np.asarray(test)
    
    # Check dimensions
    if ref_array.ndim != 2 or test_array.ndim != 2:
        raise ValueError("Input data must be 2-dimensional (samples x taxa)")
    
    if ref_array.shape[1] != test_array.shape[1]:
        raise ValueError("Reference and test datasets must have the same number of taxa")
    
    if ref_array.shape[0] == 0 or test_array.shape[0] == 0:
        raise ValueError("Input datasets must contain at least one sample")
    
    return ref_array, test_array

def load_example_data():
    """
    Load or generate example data for demonstration.
    
    Returns
    -------
    tuple
        (ref_data, test_data) as numpy arrays
    """
    # Generate synthetic microbial community data
    np.random.seed(42)
    
    # Reference dataset: 50 samples, 100 taxa
    n_ref, n_test, n_taxa = 50, 20, 100
    
    # Create reference data with some structure
    ref_data = np.random.gamma(shape=2.0, scale=2.0, size=(n_ref, n_taxa))
    # Add sparsity (many zeros) like real microbiome data
    ref_data[ref_data < np.percentile(ref_data, 40)] = 0
    
    # Create test data with slight shift
    test_data = np.random.gamma(shape=2.5, scale=1.8, size=(n_test, n_taxa))
    test_data[test_data < np.percentile(test_data, 40)] = 0
    
    return ref_data, test_data

def persistence_analysis(data, thresholds=None):
    """
    Analyze taxon persistence across different thresholds.
    
    Parameters
    ----------
    data : array-like, shape (n_samples, n_taxa)
        Microbial community data
    thresholds : array-like, optional
        Persistence thresholds to evaluate
        
    Returns
    -------
    dict
        Persistence analysis results
    """
    if thresholds is None:
        thresholds = np.linspace(0.1, 0.9, 9)
    
    data = np.asarray(data)
    persistence = np.mean(data > 0, axis=0)
    
    results = {
        'persistence': persistence,
        'taxa_persisted': [np.sum(persistence >= t) for t in thresholds],
        'thresholds': thresholds
    }
    
    return results