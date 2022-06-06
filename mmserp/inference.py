import contextlib
import time
import sys, traceback
import pkg_resources

import numpy as np
import trimesh
from scipy.spatial.distance import pdist, cdist, squareform

import pickle

from mmserp.surrogates import erp_predict
from mmserp.utils import printProgBar

STAN_MODEL_PKL_TOPHAT = pkg_resources.resource_filename(__name__, 'data/stanmodel_tophat.pkl')
STAN_MODEL_PKL_GAUSSIAN = pkg_resources.resource_filename(__name__, 'data/stanmodel_gaussian.pkl')

DATA_PATH = pkg_resources.resource_filename(__name__, 'data/surrogate_coefficients.npz')
surrogate_data = np.load(DATA_PATH)
coeffs_ERPS1 = surrogate_data["coeffs_ERPS1"]
coeffs_ERPS2 = surrogate_data["coeffs_ERPS2"]


def SD_RBF(Q, rho, alpha, scale_factor):
    """Spectral Density for RBF.
    
    Depends on dimension D, which appears as D/2, so neglected here.
    Using input Q = w^2 i.e. accepting w^2 as input, rather than w

    """
    
    SD = alpha**2 * scale_factor * (2.0 * np.pi * rho**2) * np.exp( -2 * np.pi**2 * rho**2 * Q );

    return SD


def make_s2list(s2resolution):
    """Create S2 list"""
    
    ERPmin = 160
    ERPmax = 280
    S2list_erps1 = np.arange(ERPmin, ERPmax + 1, s2resolution).astype(int)
    S2list_erps2 = np.arange(ERPmin, ERPmax + 1, s2resolution).astype(int)
    
    return S2list_erps1, S2list_erps2


def design_X(X, Tri, num, keep, NUM_designs = 100000):
    """
    Design a set of spread out vertices for assigning observations.
    
    X: coordinates
    Tri: elements (zero index)
    num: how many points
    keep: how many designs to keep
    NUM_designs: how many designs to generate and consider
    
    Summary:
    1) remove vertices too near an edge (unrealistic locations for observations)
    2) create lots of random designs and calculate minimum (smallest) distance within each
    3) keep 'keep' number of desings
    
    So we will keep 'keep' designs of 'num' points out of 'NUM_designs' possible designs.
    
    """

    # remove vertices at the edge of the mesh
    mesh = trimesh.Trimesh(vertices = X, faces = Tri, process = False)
    edges = mesh.edges_unique
    unique, counts = np.unique(mesh.faces_unique_edges, return_counts = True)
    args = np.unique(edges[unique[counts == 1]]) # this gives me the vertices in the edge
    args = np.argwhere( cdist(X, X[args]) < 6000 ) # all vertices near the edge
    
    # remove 'args' (edge vertices) from the vertices list
    verts_nonedge = np.setdiff1d(np.arange(X.shape[0]), args)
    X_nonedge = X[verts_nonedge]
    
    #print("Creating maximin design of random vertex locations...")
    dist_list = np.empty([NUM_designs])
    vert_list = np.empty([NUM_designs, num], dtype = np.int32)
    arange = np.arange(0, verts_nonedge.shape[0])
    
    for i in range(NUM_designs):
        
        # random observation vertices
        #vertIdx = np.unique( np.random.uniform(0, verts_nonedge.shape[0], num).astype(np.int32) )
        vertIdx = np.random.choice(arange, num)
        #print(vertIdx.shape)
        
        # record design
        dist_list[i] = pdist(X_nonedge[vertIdx]).min() # record minimal (Euclidean) distance
        vert_list[i] = vertIdx
    
    # order of best designs
    best = np.argsort(dist_list)[-keep:]
    tmp = vert_list[best] # this is an index into the vertices far from the edge, not an index into all vertices
    res = verts_nonedge[tmp] # should be an index into the original vertices
    
    #print("res:", res)

    return res



def stan_fit(vert_idx, erps1, erps2, S2list_erps1, S2list_erps2, deltaS2, newQ, PHI, numX, prnt = False, top_hat = True, ITER = 2000, control = {}, verbose = True, optimize = False):
    """Function for the whole stan fitting process."""

    if top_hat:
        sm = pickle.load(open(STAN_MODEL_PKL_TOPHAT, 'rb'))
    else:
        sm = pickle.load(open(STAN_MODEL_PKL_GAUSSIAN, 'rb'))
    
    vertex = []
    a_list, b_list = [],[] # collect centre of brackets, for plotting

    # what are our ERP observations for giving to STAN?
    for vv, v in enumerate(vert_idx):
        
        try:
            
            if top_hat: # new top_hat llh method
                
                a = S2list_erps1[(S2list_erps1 <= erps1[v])][-1]
                ERP_obs_a = (a)
                #print("ERP(S1) bracket", a, a + deltaS2, "for {:3.1f}".format(erps1[v]), "so ERP_obs:", ERP_obs)
                #y_erps1[vv] = ERP_obs

                b = S2list_erps2[(S2list_erps2 <= erps2[v])][-1]
                ERP_obs_b = (b)
            
            else: # original Normal llh method
                       
                a = S2list_erps1[(S2list_erps1 <= erps1[v])][-1]
                ERP_obs_a = (a + deltaS2/2.0)
                #print("ERP(S1) bracket", a, a + deltaS2, "for {:3.1f}".format(erps1[v]), "so ERP_obs:", ERP_obs)
                #y_erps1[vv] = ERP_obs

                b = S2list_erps2[(S2list_erps2 <= erps2[v])][-1]
                ERP_obs_b = (b + deltaS2/2.0)
                #print("ERP(S2) bracket", b, b + deltaS2, "for {:3.1f}".format(erps2[v]), "so ERP_obs:", ERP_obs)
                #y_erps2[vv] = ERP_obs

            a_list.append(ERP_obs_a)
            b_list.append(ERP_obs_b)
            vertex.append(v)

        except:
            pass
              
    y_erps1, y_erps2 = np.array(a_list), np.array(b_list)
    #np.zeros(len(vert_idx)), np.zeros(len(vert_idx))
        
    #print(y_erps1)
    #print(y_erps2)


    # run ERP inference with stan

    # reduce basis size while testing
    #num = NUMBASIS
    #PHI = V[:, 0:num]
    #newQ = Q[0:num]

    # NOTE: new change: set number of sub-intervals to deltaS2

    if top_hat: # new top_hat llh
        ERP_dat = {'N': len(vertex), # account for 2x observation (two ERP) in the stan code
                  'y': np.hstack([y_erps1, y_erps2]),
                  'deltaS2': deltaS2, # 95% probability that true ERP is in bracket
                  'num': int(deltaS2), # sub-intervals for top-hat approximation
                  'M': PHI.shape[1],
                  'eigen': PHI[vertex],
                  'Q': newQ,
                  'S': coeffs_ERPS1.shape[0],
                  'surrogate1': coeffs_ERPS1,
                  'surrogate2': coeffs_ERPS2}
        
    else: # original Normal llh
        ERP_dat = {'N': len(vertex), # account for 2x observation (two ERP) in the stan code
                  'y': np.hstack([y_erps1, y_erps2]),
                  'sigma': deltaS2/4.0, # 95% probability that true ERP is in bracket
                  'M': PHI.shape[1],
                  'eigen': PHI[vertex],
                  'Q': newQ,
                  'S': coeffs_ERPS1.shape[0],
                  'surrogate1': coeffs_ERPS1,
                  'surrogate2': coeffs_ERPS2}


    if not optimize:  # sampling
        fit = sm.sampling(data=ERP_dat, iter=ITER, chains=8, control = control, thin = int(2*ITER/100), verbose = verbose)
        if prnt: print(fit)
        samples = fit.extract()
    else:
        # need a while loop to restart if failed
        while True:
            try:
                # NOTE: I should choose iter for optimization carefully
                samples = sm.optimizing(data=ERP_dat, verbose = verbose, iter=ITER)
                if prnt: print(samples)
                break
            except RuntimeError as e:
                pass
                
    
    #samples = fit.extract()
    #print("sample size:", samples['beta1'].shape)
    
    with open('samples.pkl', 'wb') as f:
        pickle.dump(samples, f)
    
    beta1, beta2 = samples['beta1'], samples['beta2']
    mean1, mean2 = samples['mean1'], samples['mean2']
    rho1, rho2 = samples['rho1'], samples['rho2']
    alpha1, alpha2 = samples['alpha1'], samples['alpha2']

    # if optimization mode, need to reshape these
    if optimize:
        beta1, beta2 = beta1.reshape(1, -1), beta2.reshape(1, -1)
        mean1, mean2 = np.array([mean1]), np.array([mean2])
        rho1, rho2 = np.array([rho1]), np.array([rho2])
        alpha1, alpha2 = np.array([alpha1]), np.array([alpha2])
    
    #print(beta1, mean1)
    #import IPython
    #IPython.embed()
    
    # for ease in validation loop, I'll calculate the posterior mean here...
    # ----------------------------------------------------------------------
    
    newPhi = PHI[0:numX].T.copy()
    scale_factor = 1.0 / np.abs(PHI[0,0])
    
    f1 = np.zeros(numX)
    f2 = np.zeros_like(f1)

    count = 1
    for i in range(0, beta1.shape[0]): # use all samples remaining (after discard warm-up and thinning)

        newf1 = mean1[i] + ( beta1[i] * np.sqrt(SD_RBF(newQ, rho1[i], alpha1[i], scale_factor)) ).dot(newPhi)
        newf2 = mean2[i] + ( beta2[i] * np.sqrt(SD_RBF(newQ, rho2[i], alpha2[i], scale_factor)) ).dot(newPhi)

        f1 = f1 + (newf1 - f1) / count
        f2 = f2 + (newf2 - f2) / count

        count += 1
 
    return f1, f2, samples
    

def samples_stats(samples, newQ, PHI, numX, A):
    """Calculate statistics of fields from the samples"""

    beta1, beta2 = samples['beta1'], samples['beta2']
    mean1, mean2 = samples['mean1'], samples['mean2']
    rho1, rho2 = samples['rho1'], samples['rho2']
    alpha1, alpha2 = samples['alpha1'], samples['alpha2'] 

    # assume optimization mode based on ndim
    optimize = True if beta1.ndim == 1 else False
    # if optimization mode, need to reshape these
    if optimize:
        beta1, beta2 = beta1.reshape(1, -1), beta2.reshape(1, -1)
        mean1, mean2 = np.array([mean1]), np.array([mean2])
        rho1, rho2 = np.array([rho1]), np.array([rho2])
        alpha1, alpha2 = np.array([alpha1]), np.array([alpha2])

    print("Calculating stats of fields from {:d} samples".format(beta1.shape[0]))

    # Below, I have done this sampling to get the mean and stdev of the posterior at each point:
    # stan samples -> parameter fields samples -> erp samples -> erp mean and standard deviation.

    # for parameter sample stats
    f1 = np.zeros(numX)
    f2 = np.zeros_like(f1)
    s1 = np.zeros_like(f1)
    s2 = np.zeros_like(f1)

    # for erp samples stats
    ef1 = np.zeros_like(f1)
    ef2 = np.zeros_like(f1)
    es1 = np.zeros_like(f1)
    es2 = np.zeros_like(f1)
    

    newPhi = PHI[0:numX].T.copy()
    scale_factor = 1.0 / np.abs(PHI[0,0])
    count = 1
    #numSamples = 100
    #for i in range(beta1.shape[0] - numSamples, beta1.shape[0]):
    prefix = "progress"
    printProgBar(0, beta1.shape[0], prefix = prefix, suffix = '')
    for i in range(0, beta1.shape[0]): # use all samples remaining (after discard warm-up and thinning)
    #numSamples = beta1.shape[0]
    #for i in range(numSamples):

        # need to 'sanitize' the parameters again, because we would chuck them away if violating bounds
       
        # if I do this, do I want to calculate mean and variance of sanitized samples, or original samples?
        
        newf1 = mean1[i] + ( beta1[i] * np.sqrt(SD_RBF(newQ, rho1[i], alpha1[i], scale_factor)) ).dot(newPhi)
        newf2 = mean2[i] + ( beta2[i] * np.sqrt(SD_RBF(newQ, rho2[i], alpha2[i], scale_factor)) ).dot(newPhi)
        
        new_ef1, new_ef2, newf1, newf2 = erp_predict(newf1, newf2, A) # note: sanitizes
        
        
        tmp = f1 + (newf1 - f1) / count # tmp is updated mean
        s1 = s1 + (newf1 - f1)*(newf1 - tmp) # calculate updated var, needs old mean f1 and updated mean tmp
        f1 = tmp # save updated mean
        
        tmp = f2 + (newf2 - f2) / count
        s2 = s2 + (newf2 - f2)*(newf2 - tmp)
        f2 = tmp
        
        tmp = ef1 + (new_ef1 - ef1) / count
        es1 = es1 + (new_ef1 - ef1)*(new_ef1 - tmp)
        ef1 = tmp
        
        tmp = ef2 + (new_ef2 - ef2) / count
        es2 = es2 + (new_ef2 - ef2)*(new_ef2 - tmp)
        ef2 = tmp
        
        count += 1
        printProgBar(count, beta1.shape[0], prefix = prefix, suffix = '')
        
    # calculate stdev
    s1 = np.sqrt(s1 / (count - 1))
    s2 = np.sqrt(s2 / (count - 1))

    es1 = np.sqrt(es1 / (count - 1))
    es2 = np.sqrt(es2 / (count - 1))

    return f1, f2, s1, s2, ef1, ef2, es1, es2

