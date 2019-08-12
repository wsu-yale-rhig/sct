## A few things to note

1) In MultiplicityModel::multiplicity(double npart, double ncoll), when applying the efficiency, it applies 1/2 the efficiency for twice the ideal multiplicity. This has some effect on both tails of the distribution that seem unphysical - it slightly increases the high multiplity tail slope, and enhances the entries in the first bin (refmult=1). The reason this is done is to be consistent with StGlauber.
2) The NBDFit class requests the multiplicity with a non-integer npart, which is unphysical. This is done also to stay consistent with StGlauber.
3) For the NBD parameter optimization, it assumes no error on the glauber - but its a randomly sampled histogram, not a function, so we should be using a double-sided chi2 like what ROOT uses in its histogram chi2 test
4) in StFastGlauberMcMaker, it checks the address of an array (mIsDeformed) expecting a boolean, always returning true