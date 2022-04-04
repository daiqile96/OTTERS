#!/usr/bin/env python

"""
Markov Chain Monte Carlo (MCMC) sampler for polygenic prediction with continuous shrinkage (CS) priors.

"""

import scipy as sp
from scipy import linalg 
from scipy import random
import gigrnd

def mcmc(a, b, phi, sst_dict, n, ld_blk, blk_size, n_iter, n_burnin, thin, seed):
    print('... MCMC ...')

    # seed
    if seed != None:
        random.seed(seed)

    # derived stats
    beta_mrg = sp.array(sst_dict['BETA'], ndmin=2).T
    n_pst = (n_iter-n_burnin)/thin
    p = len(sst_dict['BETA'])

    # initialization
    beta = sp.zeros((p,1))
    psi = sp.ones((p,1))
    sigma = 1.0
    if phi == None:
        phi = 1.0; phi_updt = True
    else:
        phi_updt = False

    beta_est = sp.zeros((p,1))
    psi_est = sp.zeros((p,1))
    sigma_est = 0.0
    phi_est = 0.0

    # MCMC
    for itr in range(1,n_iter+1):
        if itr % 100 == 0:
            print('--- iter-' + str(itr) + ' ---')

        quad = 0.0
        idx_blk = range(0,blk_size)
        dinvt = ld_blk+sp.diag(1.0/psi[idx_blk].T[0])
        dinvt_chol = linalg.cholesky(dinvt)

        beta_tmp = linalg.solve_triangular(dinvt_chol, beta_mrg[idx_blk], trans='T') + sp.sqrt(sigma/n)*random.randn(len(idx_blk),1)
        beta[idx_blk] = linalg.solve_triangular(dinvt_chol, beta_tmp, trans='N')
        quad += sp.dot(sp.dot(beta[idx_blk].T, dinvt), beta[idx_blk])

        err = max(n/2.0*(1.0-2.0*sum(beta*beta_mrg)+quad), n/2.0*sum(beta**2/psi))
        sigma = 1.0/random.gamma((n+p)/2.0, 1.0/err)

        delta = random.gamma(a+b, 1.0/(psi+phi))

        for jj in range(p):
            psi[jj] = gigrnd.gigrnd(a-0.5, 2.0*delta[jj], n*beta[jj]**2/sigma)
        psi[psi>1] = 1.0

        if phi_updt == True:
            w = random.gamma(1.0, 1.0/(phi+1.0))
            phi = random.gamma(p*b+0.5, 1.0/(sum(delta)+w))

        # posterior
        if (itr>n_burnin) and (itr % thin == 0):
            beta_est = beta_est + beta/n_pst
            psi_est = psi_est + psi/n_pst
            sigma_est = sigma_est + sigma/n_pst
            phi_est = phi_est + phi/n_pst

    # print estimated phi
    if phi_updt == True:
        print('... Estimated global shrinkage parameter: %1.2e ...' % phi_est )

    print('... Done ...')

    return beta_est
