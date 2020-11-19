from __future__ import print_function
from orphics import mpi
import numpy as np
import os,sys
from classy import Class
import camb
from camb import model

"""
Unified framework for supporting:
1. CMB lensing (CLASS or CAMB)
2. CMB power spectra (CLASS or CAMB)
3. BAO (CLASS or CAMB)
4. hmvec short wavelength Pge, Pgg
5. long wavelength Pgv, Pgg, Pvv

"""

# def bandpower_gaussian_covariance


# def save_symm_derivs(,fiducials,step_fracs,output_root):
#     """
#     Save symmetric derivatives for the given observable.

#     """

#     fs = fiducials
#     vparams = sorted(step_fracs.keys())
#     assert len(vparams)==len(set(vparams))

#     nparams = len(vparams)

#     comm,rank,my_tasks = mpi.distribute(nparams)

#     for task in my_tasks:
#         param = vparams[task]

#         step = step_sracs[param] * fs[param]

#         fparams = dict(self.fs)
#         fparams[param] = fparams[param] + step
#         forward = self.get_obs(observables,fparams)

#         fparams = dict(self.fs)
#         fparams[param] = fparams[param] - step
#         backward = self.get_obs(observables,fparams)

#         for obs in observables:
#             fname = output_root + f"_{obs}_{param}_deriv.txt"
#             d = (forward[obs] - backward[obs])/2./step

        

def map_params(params,engine='camb'):
    h = params['H0'] / 100.
    if 'omm' in params.keys():
        assert 'omch2' not in params.keys()
        

    if engine=='camb':
        return params
    elif engine=='class':
        return _camb_to_class(params)


def set_defaults(params):
    ds = {
        'omch2': 0.1203058
        ,'ombh2': 0.02219218
        ,'H0': 67.02393
        ,'ns': 0.9625356
        ,'As': 2.15086031154146e-9
        ,'mnu': 0.06
        ,'w0': -1.0
        ,'wa': 0.0
        ,'tau':0.06574325
        ,'nnu':3.046
        ,'ok':0
        ,'r':0
        ,'cs2':1.0
        ,'thetastar': None
        ,'ctheta': None
    }
    for key in ds.keys():
        if key not in params.keys(): params[key] = ds[key]
    return params

    
    

def set_camb_pars(params=None,de='fluid'):
    """
    We only allow an As,ns,omch2,ombh2,H0,tau,mnu,w0,wa,omk,r,ctheta,nnu,
    cs2,thetastar parametrization. (ctheta=cosmomc_theta)
    Any other parameterization has to be obtained by transforming the Fisher matrix.
    """
    if params is None: params = dict({})
    params = set_defaults(params)
    pars = camb.CAMBparams()
    #This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
    pars.set_cosmology(H0=params['H0'], ombh2=params['ombh2'], 
                       omch2=params['omch2'], mnu=params['mnu'], 
                       omk=params['ok'], tau=params['tau'],nnu=params['nnu'],
                       cosmomc_theta=params['ctheta'],thetastar=params['thetastar'])
    pars.InitPower.set_params(As=params['As'], ns=params['ns'], r=params['r'])
    pars.set_dark_energy(w=params['w0'], wa=params['wa'], cs2=params['cs2'],dark_energy_model=de)
    return pars


def get_s8(zs=[0.],params=None,nonlinear=False,kmax=5.2,**kwargs):
    zs = np.asarray(zs)
    zdiffs = np.diff(zs)
    if np.all(zdiffs>0):
        rev = True
    elif np.all(zdiffs<0):
        rev = False
    else:
        raise ValueError
    pars = set_camb_pars(params=params,**kwargs)
    pars.set_matter_power(redshifts=zs,kmax=kmax,nonlinear=nonlinear,silent=True)
    if nonlinear: 
        pars.NonLinear = model.NonLinear_both
    else:
        pars.NonLinear = model.NonLinear_none
    results = camb.get_results(pars)
    s8 = np.array(results.get_sigma8())
    return s8[::-1] if rev else s8

def get_bao(params,zs,engine='camb'):
    params = map_params(params,engine=engine)
    if engine=='camb':
        pars = camb.CAMBparams()
        #This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
        pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
        pars.InitPower.set_params(As=2e-9, ns=0.965, r=0)
        pars.set_for_lmax(2500, lens_potential_accuracy=0)
    elif engine=='class':
        cosmo = Class()
        params['output'] = ''
        cosmo.set(params)
        cosmo.compute()
        Hzs = np.array([cosmo.Hubble(z) for z in zs])
        D_As = np.array([cosmo.angular_distance(z) for z in zs])
        D_Vs = ((1+zs)**2 * D_As**2 * zs/Hzs)**(1/3.)
        retval = cosmo.rs_drag()/D_Vs
        cosmo.struct_cleanup()
        cosmo.empty()
    return retval


def get_cls(params,lmax=3000,engine='camb'):
    params = map_params(params,engine=engine)
    if engine=='camb':
        pass
    elif engine=='class':
        cosmo = Class()
        params['output'] = 'lCl,tCl,pCl'
        params['lensing'] = 'yes'
        params['non linear'] = 'hmcode'
        params['l_max_scalars'] = lmax
        cosmo.set(params)
        cosmo.compute()
    return retval


def get_trans_deriv(iparam,oparam,fiducials,deriv_path):
    pass
    

def reparameterize(Fmat,iparams,oparams,fiducials,deriv_path):
    """
    Re-parameterize a Fisher matrix.
    iparams must only contain CAMB primitives
    like As,ns,omch2,ombh2,H0,tau,mnu,w0,wa,omk,r.
    oparams can contain those as well as
    s8
    om_m
    om_b
    Symbol

    where Symbol can be a sympy.Symbol expression,
    e.g. Symbol(s8) * (Symbol(om_m)/0.3)**0.25
    """

    onum = len(oparams)
    inum = len(iparams)
    
    M = np.zeros((inum,onum))

    for i in range(inum):
        for j in range(onum):
            M[i,j] = get_trans_deriv(iparams[i],oparams[j],fiducials,deriv_path)

    return np.dot(np.dot(M,Fmat),M.T)


#def s8(params):
    
