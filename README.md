# SCR-Multisession
Basic multisession SCR MCMC sampler in nimble using RJMCMC instead of data augmentation. Poisson observation model.

These models use count prior data augmentation for group size only. Could be used for number of groups as well. 
https://github.com/benaug/SCR-Count-Prior-Data-Augmentation

See Royle and Converse (2014) for another way to fit multisession SCR models while retaining the Poisson assumption.
https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12135