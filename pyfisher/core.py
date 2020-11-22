import camb
import numpy as np

def getHubbleCosmology(theta,params,precision = 0.00001,H0init = 70.,H0max = 100.,H0min = 40.):
    '''
    Finds H0 for given theta and other cosmo params{} using a bisection search
    '''
    
    H0 = H0init
    err = precision*2.
    pars = camb.CAMBparams()
    pars.set_accuracy(AccuracyBoost=2.0, lSampleBoost=2.0, lAccuracyBoost=2.0)
    #pars.set_accuracy(AccuracyBoost=3.0, lSampleBoost=3.0, lAccuracyBoost=3.0)
    pars.set_dark_energy(w=params['w'],wa=params['wa'],dark_energy_model='ppf')
    pars.InitPower.set_params(As=params['As'],ns=params['ns'],r=params['r'])
    #no need for this: pars.set_for_lmax(lmax=int(params['lmax']), lens_potential_accuracy=1, max_eta_k=2*params['lmax'])

    print "Searching for H0 through bisection..."
    j = 0
    while np.abs(err)>precision:
        pars.set_cosmology(H0=H0, ombh2=params['ombh2'], omch2=params['omch2'], tau=params['tau'],mnu=params['mnu'],nnu=params['nnu'],omk=params['omk'])
        results = camb.get_results(pars)
        HundredTheta = 100.*results.cosmomc_theta()

        err = HundredTheta-theta

        if err<0.:
            H0min = H0
            H0 = (H0+H0max)/2.
        elif err>0.:
            H0max = H0
            H0 = (H0+H0min)/2.

        j+=1
        

    print "Found H0 = ", H0, " in ", j , " iterations."
    return H0
