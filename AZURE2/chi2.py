import os
import array

import numpy as np

from pyazr import azure2
from lmfit import Parameters, Minimizer

os.environ["OMP_NUM_THREADS"] = "1"

# Nuisance parameter map : { index : (value, error) } (in AZURE2 order)
nuisances = {  }

# We read the .azr file
azr = azure2('7Be.azr', nprocs=1)

# We get the initial values from AZURE2
theta0 = azr.params
ntheta = len(theta0)

# We'll read the data from the output file since it's already in the center-of-mass frame
y = azr.cross
yerr = azr.cross_err

# Calculate number of data points
ndata = 0
for i in range(len(y)):
    ndata += len(y[i])

# Now add the normalization parameters to theta
norms = {  }
for i, data in enumerate(y):
    theta0 = np.concat( (theta0, [1.0]) )
    norms[i] = (1.0, 0)

# Callback function to print the chi2 at each iteration
def callback(params, iter, resid, *args, **kws):
    if( iter % 50 == 0 ): 
        print("Iteration: {:8d} it Chi2: {:15.4f}".format( iter, np.sum( np.array(resid)**2 ) ), end="\r")
    pass

# Add nuisance parameter to chi2
def nuisance( theta ):
    nu = []
    for i in nuisances: nu.append( (theta[i] - nuisances[i][0]) / nuisances[i][1] )
    return nu

# Add nuisance normalization to chi2
def normalization( theta ):
    norm = []
    for i in norms:
        if( norms[i][1] != 0 ):
            idx = ntheta + i
            norm.append( (theta[idx] - norms[i][0]) / norms[i][1] )
    return norm

# Calculated squared residuals
def least_squares( mu, theta ):
    res = []
    for i in range( len( mu ) ):
        idx = ntheta + i
        res.extend( (mu[i] - y[i] * theta[idx]) / ( yerr[i] * theta[idx] ) )
    return res

#Function to minimize
def func( theta ):
    theta = list( theta.valuesdict().values() )
    mu = azr.calculate( theta[:ntheta], proc=0 )
    fcn = least_squares( mu, theta )
    fcn.extend( normalization( theta ) )
    fcn.extend( nuisance( theta ) )
    return fcn

# Create the parameters
params = Parameters()
for i, param in enumerate(theta0):
    params.add('param_{}'.format(i), value=param, vary=True)

# Starting the minimization
minimizer = Minimizer(func, params, nan_policy='raise', iter_cb=callback)
out = minimizer.minimize(method='nelder')

# Getting the parameters
result = [param.value for _, param in out.params.items()]
covari = out.covar

#  Get the final chi2
chi2 = np.sum( pow( out.residual, 2 ) )
chi2_red = chi2 / (ndata - len(theta0))
print("Chi2: {:15.4f}".format( chi2 ))
print("Chi2 Reduced: {:15.4f}".format( chi2_red ))
print("Parameters: ")
for i, param in enumerate(result):
    try:
        print("  {:2d} : {:15.4f} +/- {:15.4f}".format( i, param, np.sqrt(covari[i][i]) ))
    except:
        print("  {:2d} : {:15.4f} +/- {:15.4f}".format( i, param, 0.0 ))

# Write the results to a file
with open("results/results_chi2.txt", "w") as f:
    f.write("Chi2: {:15.4f}\n".format( chi2 ))
    f.write("Chi2 Reduced: {:15.4f}\n".format( chi2_red ))
    f.write("Parameters:\n")
    for i, param in enumerate(result):
        try:
            f.write("  {:2d} : {:15.4f} +/- {:15.4f}\n".format( i, param, np.sqrt(covari[i][i]) ))
        except:
            f.write("  {:2d} : {:15.4f} +/- {:15.4f}\n".format( i, param, 0.0 ))

# Write a parameter file
with open("results/parameters_chi2.txt", "w") as f:
    for i, param in enumerate(result):
        f.write("{:15.4f}\n".format( param ))

# Write a covariance file
#with open("results/covariance_chi2.txt", "w") as f:
#    for i in range(len(covari)):
#        for j in range(len(covari[i])):
#            f.write("{:15.4f}".format( covari[i][j] ))
#        f.write("\n")

# Exit the program
print("Finished!")