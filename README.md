# RPD
Function and Analysis Scripts for Ratio Percentile Deviation Method:

The RPD function is intended to be used as follows: copy or download the RPD function script, and run the entire script to load the calc.rpd() function into your R environment. Then, use the function as follows:

calc.rpd(ref = XXXX, test = YYYY, ref.pers = ZZZZ)

An optional argument allows the user to also retrieve the taxon scores:

calc.rpd(ref = XXXX, test = YYYY, ref.pers = ZZZZ, taxon.scores = T)

The RPD function R script has supporting functions and a main calc.rpd function. The supporting functions are run within the calc.rpd function and do not need to be run by the user. 

A package for Python does not yet exist, but you may find success in asking Claude code to translate the function from R to Python

FAQ:
- How many samples do I need in the reference set? In the test set?
  If working with experimental data where samples are homogenous, you may be able to use as few as 8-10 samples for a reference.
  The exact threshold is unknown without more empirical data, but this seems like a reasonable extrapolation from the observation that even 20 references worked well for a very heterogeneous dataset.
  There is no limit on test samples, although the RPD function does not handle single samples in the test set due to the reliance on matrix functions.
  If needing to calculate RPD on a single sample, simply duplicate the sample in a matrix or run it with other samples. 
- How do I use RPD with experimental data?
  The most important thing is to have a representative dataset as a reference. This is likely the microbial community before perturbation, or in the control treatment.
  For each treatment, use those samples as the test dataset and the initial conditions or control treatment as the reference. You can then run statistical analyses on RPD values just like you would on diversity values. 
  It is best to use the same reference persistence threshold for all RPD calculations using the same reference if you intend to compare RPD values!
- How do I select an appropriate reference persistence threshold?
  Your analysis should ideally use a reference persistence threshold where minor changes (i.e. 0.05 either direction) do not strongly change the results.
  For deeply sequenced datasets with high diversity, the persistence threshold will likely be higher.
  I generally started at 0.7 and then compared against 0.5 and 0.85, then moved up or down depending on which of those options was better.
  You can also try the analyses shown in the supplementary materials to more rigorously choose a reference persistence threshold. 
