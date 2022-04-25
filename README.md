# SCR-Multisession
Basic multisession SCR MCMC sampler in nimble using RJMCMC instead of data augmentation. Poisson observation model.

Using reversible jump MCMC instead of data augmentation allows one to retain the Poisson assumption N[g] ~ dpois(lambda[g]) for session g.

Using data augmentation, the Poisson assumption for N[g] is replaced by N[g] ~ dbinom(psi[g],M[g]). This converges to a Poisson when M[g] -> infinity keeping N[g] fixed. Fitting the model this way would be incredibly computationally inefficient.

See Royle and Converse (2014) for another way to fit multisession SCR models while retaining the Poisson assumption.
https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12135