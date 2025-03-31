# TransPairs

The code in this repository was used to perform the transmission analysis described in the following paper:

Eldholm et al (2016) Impact of HIV co-infection on the evolution and transmission of multidrug-resistant tuberculosis. 
Elife 5:e16644
http://elifesciences.org/lookup/doi/10.7554/eLife.16644

-   The `anaPairs.m` file is a Matlab script computing the pairwise matrix of likelihoods of transmission from any host to any other.
-   The `edmunds.R` file is a R script that extracts from this matrix the most likely transmission events.
-   The `TransPairs.R` file is a newer R script equivalent to both scripts above.

The comments in the files provide more information about how to use them. 
The Methods section of the paper provides an explanation of the methodology implemented here.
