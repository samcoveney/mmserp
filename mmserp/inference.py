import contextlib
import time
import sys, traceback
import pkg_resources

STAN_MODEL_PKL = pkg_resources.resource_filename(__name__, 'data/stanmodel_tophat.pkl')

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



def stan_fit(vert_idx, erps1, erps2, S2list_erps1, S2list_erps2, deltaS2, prnt = False, control = {}):
    """Function for the whole stan fitting process."""

    # reload this to rededine the object just in case
    sm = pickle.load(open(STAN_MODEL_PKL, 'rb'))
    
    # good for testing
    ITER = 2000

    vertex = []
    a_list, b_list = [],[] # collect centre of brackets, for plotting

    # what are our ERP observations for giving to STAN?
    for vv, v in enumerate(vert_idx):
        
        try:
            
            if True: # new top_hat llh method
                
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


    if True: # new top_hat llh
        ERP_dat = {'N': len(vertex), # account for 2x observation (two ERP) in the stan code
                  'y': np.hstack([y_erps1, y_erps2]),
                  'deltaS2': deltaS2, # 95% probability that true ERP is in bracket
                  'num': 10, # sub-intervals for top-hat approximation
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


    fit = sm.sampling(data=ERP_dat, iter=ITER, chains=4, control = control, thin = 10)

    if prnt: print(fit)
    
    samples = fit.extract()
    print("sample size:", samples['beta1'].shape)
    
    #with open('samples.pkl', 'wb') as f:
    #    pickle.dump(samples, f)
    
    beta1, beta2 = samples['beta1'], samples['beta2']
    mean1, mean2 = samples['mean1'], samples['mean2']
    rho1, rho2 = samples['rho1'], samples['rho2']
    alpha1, alpha2 = samples['alpha1'], samples['alpha2']

    
    # for ease in validation loop, I'll calculate the posterior mean here...
    # ----------------------------------------------------------------------
    
    # Note that we used thin = 1, and yet I've only used last numSamples
    
    # TODO: should discard the burn-in (although Stan may have done this), thin the rest, the select random runs
    
    newPhi = PHI[0:X.shape[0]].T.copy()
    
    f1 = np.zeros(X.shape[0])
    f2 = np.zeros_like(f1)

    count = 1
    numSamples = 100
    #for i in range(beta1.shape[0] - numSamples, beta1.shape[0]):
    for i in range(0, beta1.shape[0]): # use all samples remaining (after discard warm-up and thinning)

        newf1 = mean1[i] + ( beta1[i] * np.sqrt(SD_RBF(newQ, rho1[i], alpha1[i])) ).dot(newPhi)
        newf2 = mean2[i] + ( beta2[i] * np.sqrt(SD_RBF(newQ, rho2[i], alpha2[i])) ).dot(newPhi)

        f1 = f1 + (newf1 - f1) / count
        f2 = f2 + (newf2 - f2) / count

        count += 1
 
    return f1, f2, samples
    

