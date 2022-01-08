# ERPsurrogate.py
# routines to create parameter fields and corresponding ERP maps using the ERP surrogate functions

# imports
import pkg_resources
import numpy as np
from scipy.spatial.distance import pdist, cdist, squareform
from scipy.special import gamma


# load the surrogate data and create the surrogate function

DATA_PATH = pkg_resources.resource_filename(__name__, 'data/surrogate_coefficients.npz')
surrogate_data = np.load(DATA_PATH)
coeffs_ERPS1 = surrogate_data["coeffs_ERPS1"]
coeffs_ERPS2 = surrogate_data["coeffs_ERPS2"]


def spectralDensity(w, smoothness, lengthscale, amplitude):
    """Spectral Density for Matern kernel.
    
       w: sqrt(eigenvalues)
       smoothness: Matern covariance smoothness 
       lengthscale: lengthscale parameter
       amplitude: i.e. the sigma**2 multiplying the kernel

    """
    
    D = 2  # dimension
    v = smoothness  # smoothness
    rho = lengthscale # lengthscale

    alpha = (2**D * np.pi**(D/2) * gamma(v + (D/2)) * (2*v)**v) / gamma(v)
    delta = 1.0 / rho**(2*v)
    beta = 2*v / (rho**2)  + 4*(np.pi**2)*(w**2) 
    expo = -(v + (D/2))

    result = (alpha*delta) * beta**expo

    result = result * amplitude

    return result


def SD_RBF(Q, rho, alpha):
    """Spectral Density for RBF.
    
    Depends on dimension D, which appears as D/2, so neglected here.
    Using input Q = w^2 i.e. accepting w^2 as input, rather than w

    """
    
    SD = alpha**2 * (2.0 * np.pi * rho**2) * np.exp( -2 * np.pi**2 * rho**2 * Q );

    return SD


def prior_sample(Q, V, lengthscale, amplitude):
    """Generate a draw from the prior of the GP on the manifold."""

    #print("Problem: I generate fields with 5/2 smoothness, but  I fit a model with 3/2 smoothness?")
    
    SD = spectralDensity(w = np.sqrt(Q), smoothness = 5./2., lengthscale = lengthscale, amplitude = amplitude)

    coeffs = np.random.normal(loc = 0, scale = np.sqrt(SD))

    prior = coeffs.dot(V.T)

    return prior


def surrogate_ERP(tau_out, APD_max):
    """Surrogate ERP model.

       Inputs must be np.array type

    """
    
    if tau_out.ndim == 0: ONES = 1
    else: ONES = np.ones_like(tau_out)
    
    #print("HERE:", ONES.shape, tau_out.shape)
    
    ERP_basis = np.vstack([ONES, tau_out, APD_max,
                           tau_out**2, APD_max**2, tau_out*APD_max,
                           tau_out**3, APD_max**3, tau_out**2*APD_max, tau_out*APD_max**2])

    ERPS1 = np.dot(coeffs_ERPS1, ERP_basis)
    ERPS2 = np.dot(coeffs_ERPS2, ERP_basis)

    
    if tau_out.ndim == 0:
        if ERPS1 > 285: ERPS1, ERPS2 = np.nan, np.nan
    else:
        cond = (ERPS1 > 285)
        ERPS1[cond], ERPS2[cond] = np.nan, np.nan
        ERPS1, ERPS2 = np.squeeze(ERPS1), np.squeeze(ERPS2)

    return ERPS1, ERPS2

    
def sq_biharmonic(A, k, good):
    """Biharmonic distance."""
    return np.sum( (A[good] - A[k])**2 , axis = 1)
    

def erp_predict(tau_out, APD_max, A):
    """
    Predict the full ERP field using the surrogate.
    
    Unfortunately, this function is still using some global variable... never mind for now.
    """
    
    # first, replace tau_out if below the minimum
    bad_tau = (tau_out < 1)
    tau_out[bad_tau] = 1.0

    if False:
        erps1, erps2 = np.empty(X.shape[0]), np.empty(X.shape[0])
        for i in range(X.shape[0]):
            tmp1, tmp2 = surrogate_ERP(tau_out[i], APD_max[i])
            #print(i, tau_out[i], APD_max[i], "--->", tmp1, tmp2)
            erps1[i], erps2[i] = tmp1, tmp2
            
    if True: # vectorized
        erps1, erps2 = surrogate_ERP(tau_out[None,:], APD_max[None,:])
        
    #print( np.allclose(erps1, erps1_test.flatten()) )
    #print( np.allclose(erps2, erps2_test.flatten()) )

    # replace NaN values with best non-nan values
    # based on 'nearest in space' non-NaN value
    # -------------------------------------------

    # vertex index of good and bad values
    # condition will be where stuff is bad
    # additional condition to keep  165 < ERP(S1) < 275 (170 to 270 is realistic range, so adding buffer)
    
    # The discontinuity is at 285, so perhaps I can be more flexible here...
    # on the other hand, if we really believe in erps1 < 275 for example, then keep this here
    
    condition_nan = np.isnan(erps1)
    erps1[condition_nan] = 1000
    condition_2 = (erps1 > 279)
    erps1[condition_nan] = -1000
    condition_3 = (erps1 < 161)
    
    condition = condition_nan | condition_2 | condition_3
    
    bad = np.argwhere(condition).flatten()
    good = np.argwhere(~condition).flatten()
    #print("bad:", bad)


    # use an approximate over-the-manifold distance
    
    for bb, b in enumerate(bad):

        # distance weight average method
        w = 1.0 / sq_biharmonic(A, b, good)**4

        # average the parameters, and calculate the resulting ERP
        tau_out[b] = np.sum(w * tau_out[good])/np.sum(w)
        APD_max[b] = np.sum(w * APD_max[good])/np.sum(w)
        erps1[b], erps2[b] = surrogate_ERP(tau_out[b], APD_max[b])
        
    return erps1, erps2, tau_out, APD_max
    

# create ERP feature (i.e. tau_out and APD_max) maps
def scale_range(field, MIN, MAX):
    """Scale the field into the minmax range."""
    
    # field goes from zero to one
    tmp = (field - field.min()) / (field.max() - field.min())
    
    # scale into minmax range
    tmp = MIN + (MAX - MIN)*tmp
    
    return tmp
    

# create maps of ERP over the atrium
def ERP_maps(X, Tri, Q, V, lengthscale):
    """Create random draws of ERP(S1) and ERP(S2) maps"""
    
    # MINMAX of features
    tau_out_minmax = [1.0, 30.0]
    APD_max_minmax = [120.0, 250.0]

    # required for biharmonic distance calculation
    A = V[0:X.shape[0], 1:33]/Q[1:33]
    
    # draw from prior; create sample from several samples with different lengthscales
    # NOTE: samples will be scaled into minmax range, so amplitude is irrelevant
    sample_1 = prior_sample(Q, V[0:X.shape[0]], lengthscale = lengthscale, amplitude = 10.0)
    sample_2 = prior_sample(Q, V[0:X.shape[0]], lengthscale = lengthscale, amplitude = 10.0)

    tau_out = scale_range(sample_1, tau_out_minmax[0], tau_out_minmax[1])
    APD_max = scale_range(sample_2, APD_max_minmax[0], APD_max_minmax[1])
    #features_plot_for_paper(X, Tri, tau_out, title = "tau_out (true)").show()
    #features_plot_for_paper(X, Tri, APD_max, title = "APD_max (true)").show()

    # NOTE: erp_predict will also sanitize the parameters, so that they generate valid ERPs
    erps1, erps2, tau_out, APD_max = erp_predict(tau_out, APD_max, A)
    
    return erps1, erps2, tau_out, APD_max


# create ERP samples from parameter samples, also sanitizing the parameters
def make_erp_samples(tsamples, asamples, getSamples):

    if getSamples == True:

        erps1_samples_TS = np.zeros_like(tsamples)
        erps2_samples_TS = np.zeros_like(tsamples)
        tau_out_samples_TS = np.zeros_like(tsamples)
        apd_max_samples_TS = np.zeros_like(tsamples)

        for ii in range(tsamples.shape[1]):

            erps1_samples_TS[:,ii], erps2_samples_TS[:,ii], tau_out_samples_TS[:,ii], apd_max_samples_TS[:,ii] = erp_predict(tsamples[:,ii], asamples[:,ii])

    else:

        print("[WARNING]: did not calculate full posterior covariance, so have no samples...")

    return erps1_samples_TS, erps2_samples_TS, tau_out_samples_TS, apd_max_samples_TS


