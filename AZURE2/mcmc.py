import os
import emcee

import numpy as np

from scipy import stats
from pyazr import azure2
from multiprocess import Pool, current_process

# Restrict processes to one thread only
os.environ['OMP_NUM_THREADS'] = '1'

# Define the parameters prior distributions
priors = [ ]

# Minimization variables
nsteps   = 100000       # How many steps should each walker take?
nprocs   = 5            # How many Python processes do you want to allocate?

# We read the .azr file and set the external capture file to speed up the calculation
azr = azure2('7Be.azr', nprocs=nprocs)

# We get the initial values from AZURE2
theta0 = azr.params
ntheta = len(theta0)

# We'll read the data from the output file since it's already in the center-of-mass frame
y = azr.cross
yerr = azr.cross_err

# We define variables
ndim     = len(theta0) + len(y)  # How many parameters are you fitting?
nwalkers = 2 * ndim              # How many walkers do you want to use?

# Create priors for all theta0, these must be uniforms in the 90% confidence interval
for i, theta in enumerate(theta0):
    priors.append(stats.uniform(-1e12, 2e12))

# Add the priors for data that are uniform between 0 and 10 for each data in y
for i, data in enumerate(y):
    priors.append(stats.uniform(0.1, 9.9))

# Prior log probability
def lnPi( theta ):
    return np.sum([pi.logpdf(t) for (pi, t) in zip(priors, theta)])

# Log likelihood
def lnL( theta, proc=0 ):
    res = 0
    mu = azr.calculate( theta[:ntheta], proc=proc )
    for i in range( len( mu ) ):
        idx = ntheta + i
        res += -0.5 * np.sum( np.log(2 * np.pi * pow(yerr[i], 2) ) + pow((mu[i] - y[i] * theta[idx]) / (yerr[i] * theta[idx]), 2) )
    return res

# Posterior log probability
def lnP( theta ):
    try: proc = int(current_process().name.split('-')[1]) - 1 # We want to get the numbe r of the process to call the right AZURE2 port
    except: proc = 0 
    lnpi = lnPi( theta )
    if not np.isfinite( lnpi ): return -np.inf
    return lnL( theta, proc=proc ) + lnpi

# Read best fit from the output file
best = np.loadtxt('results/parameters_chi2.txt')

# Prepare initial walker positions
p0 = np.zeros( (nwalkers, ndim) )
for i in range(nwalkers):
    for j in range(ndim):
        if( j < ntheta ): 
            if( best[j] > 0 ): p0[i, j] = np.random.normal(best[j], 1e-5 * best[j])
            else: p0[i, j] = np.random.normal(best[j], -1e-5 * best[j])
        else: p0[i, j] = np.random.normal(best[j], 1e-5 * best[j])

# Prepare the file to write the chains
backend = emcee.backends.HDFBackend('results/samples.h5') 

try:
    n_iterations_done = backend.iteration
    print(f"Restarting from iteration {n_iterations_done}")
    with Pool(processes=nprocs) as pool:
        sampler = emcee.EnsembleSampler( nwalkers, ndim, lnP, pool=pool, backend=backend ) 
        state = sampler.run_mcmc( None, nsteps, progress=True, tune=True )

except:
    print("Starting from scratch")
    with Pool(processes=nprocs) as pool:
        sampler = emcee.EnsembleSampler( nwalkers, ndim, lnP, pool=pool, backend=backend ) 
        state = sampler.run_mcmc( p0, nsteps, progress=True, tune=True )