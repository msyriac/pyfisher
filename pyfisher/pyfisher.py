from __future__ import print_function
from orphics import mpi,cosmology,stats
import numpy as np
import os,sys
from classy import Class
import camb
from camb import model
import sympy
import warnings

"""
Unified framework for supporting:
1. CMB lensing (CLASS or CAMB)
2. CMB power spectra (CLASS or CAMB)
3. BAO (CLASS or CAMB)
4. hmvec short wavelength Pge, Pgg
5. long wavelength Pgv, Pgg, Pvv

"""

def get_binner(bin_edges,interpolate):
    if interpolate:
        cents = (bin_edges[1:] + bin_edges[:-1])/2.
        bin = lambda x: x(cents)
    else:
        binner = stats.bin1D(bin_edges)
        ells = np.arange(bin_edges[0],bin_edges[-1]+1,1)
        bin = lambda x: binner.bin(ells,x(ells))[1]
        cents = binner.cents
    return cents, bin

def gaussian_band_covariance(bin_edges,specs,cls_dict,nls_dict,interpolate=False):

    cents,bin = get_binner(bin_edges,interpolate)
    delta_ell = np.diff(bin_edges)

    def _symmz(cdict,ab):
        try:
            return bin(cdict[ab])
        except KeyError:
            try:
                return bin(cdict[ab[::-1]])
            except KeyError:
                return cents*0
            
    
    ncomps = len(specs)
    nbins = len(bin_edges) - 1
    cov = np.zeros((nbins,ncomps,ncomps))
    for i in range(ncomps):
        for j in range(i,ncomps):
            spec1 = specs[i]
            spec2 = specs[j]

            a,b = spec1
            g,d = spec2

            ag = a+g
            bd = b+d
            ad = a+d
            bg = b+g

            cl_ag = _symmz(cls_dict,ag)
            nl_ag = _symmz(nls_dict,ag)
            cl_bd = _symmz(cls_dict,bd)
            nl_bd = _symmz(nls_dict,bd)
            cl_ad = _symmz(cls_dict,ad)
            nl_ad = _symmz(nls_dict,ad)
            cl_bg = _symmz(cls_dict,bg)
            nl_bg = _symmz(nls_dict,bg)

            cov[:,i,j] = ((cl_ag+nl_ag)*(cl_bd+nl_bd)+(cl_ad+nl_ad)*(cl_bg+nl_bg))/(2*cents+1)/delta_ell
            if i!=j: cov[:,i,j] = cov[:,j,i].copy()
    return cov


def band_fisher(param_list,bin_edges,specs,cls_dict,nls_dict,derivs_dict,interpolate=True):

    cents,bin = get_binner(bin_edges,interpolate)
    nbins = len(bin_edges) - 1
    ncomps = len(specs)

    cov = gaussian_band_covariance(bin_edges,specs,cls_dict,nls_dict,interpolate=interpolate)
    cinv = np.linalg.inv(cov)

    nparams = len(param_list)
    Fisher = np.zeros((nparams,nparams))
    for i in range(nparams):
        for j in range(i,nparams):

            param1 = param_list[i]
            param2 = param_list[j]
            dcls1 = np.zeros((nbins,ncomps))
            dcls2 = np.zeros((nbins,ncomps))
            for k,spec in enumerate(specs):
                dcls1[:,k] = bin(derivs_dict[param1][spec])
                dcls2[:,k] = bin(derivs_dict[param2][spec])

            Fisher[i,j] = np.einsum('ik,ik->',np.einsum('ij,ijk->ik',dcls1,cinv),dcls2)
            if i!=j: Fisher[j,i] = Fisher[i,j]

    return stats.FisherMatrix(Fisher,param_list)




def get_jobs(param_file,exclude=None):
    param_dat = np.genfromtxt(param_file,dtype=None,encoding='utf-8',delimiter=',')
    jobs = []
    jobs.append((None,None,'f'))
    fids = {}
    for p in param_dat:
        param = p[0]
        if exclude is not None:
            if param in exclude:
                print(f"Skipping {param}")
                continue
        fid = p[1]
        fids[param] = fid
        pstr = str(p[2]).strip()
        if pstr[-1]=='%': 
            step = float(pstr[:-1])*np.abs(fid)/100.
        else:
            step = float(pstr)
        assert step>0
        jobs.append((param,fid+step,'u'))
        jobs.append((param,fid-step,'d'))
    return jobs,fids

def _camb_to_class(params):
    if params['thetastar'] is not None:
        params['100*theta_s'] = params.pop('thetastar')*100.
    else:
        params.pop('thetastar')
    if params['ctheta'] is not None:
        params['100*theta_s'] = params.pop('ctheta')*100.
        warnings.warn('Replacing CLASS 100*theta_s with cosmomc_theta')
    else:
        params.pop('ctheta')
    if params['H0'] is None:
        params.pop('H0')
    params['A_s'] = params.pop('As')
    params['n_s'] = params.pop('ns')
    params['Omega_cdm'] = params.pop('omch2')
    params['omega_b'] = params.pop('ombh2')
    params['tau_reio'] = params.pop('tau')
    params['Omega_k'] = params.pop('ok')
    params['N_ur'] = params.pop('nnu')
    #params['m_ncdm'] = ','.join([str(params.pop('mnu')/3.)]*3)
    params['N_ncdm'] = 1
    params['m_ncdm'] = params.pop('mnu')
    
    # params['use_ppf'] = 'no'
    # params['fluid_equation_of_state'] = 'CLP'
    # params['w0_fld'] = params.pop('w0')
    # params['wa_fld'] = params.pop('wa')
    # params['cs2_fld'] = params.pop('cs2')

    # params.pop('mnu')
    # DARK ENERGY AND R NOT SUPPORTED
    params.pop('w0')
    params.pop('wa')
    params.pop('cs2')
    params.pop('r')
    return params


def map_params(params,engine='camb'):
    checksum = 0
    try:
        assert params['H0']
        checksum = checksum + 1
    except:
        pass
    try:
        assert params['thetstar']
        checksum = checksum + 1
    except:
        pass
    try:
        assert params['ctheta']
        checksum = checksum + 1
    except:
        pass
    assert checksum <= 1
    if params is None: params = dict({})
    params = set_defaults(params)
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

    
    

def set_camb_pars(params=None,de='ppf'):
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
    if params['r']>0: pars.WantTensors = True
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

def get_bao_rs_dV(zs,params=None,engine='camb',de='ppf'):
    #FIXME: camb and class only agree at 3% level!!!
    params = map_params(params,engine=engine)
    if engine=='camb':
        pars = set_camb_pars(params=params,de=de)
        results = camb.get_results(pars)
        retval = results.get_BAO(zs,pars)[:,0]
    elif engine=='class':
        zs = np.asarray(zs)
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


def get_bao_fisher_rs_dV_diagonal(param_list,deriv_theory,fiducial_theory,sigma_percents):
    nparams = len(param_list)
    Fisher = np.zeros((nparams,nparams))
    sigmas = sigma_percents*fiducial_theory/100.
    for i in range(nparams):
        for j in range(i,nparams):
            param1 = param_list[i]
            param2 = param_list[j]
            rsdV1 = deriv_theory[param1]
            rsdV2 = deriv_theory[param2]
            Fz = (rsdV1*rsdV2/sigmas**2.).sum()
            Fisher[i,j] = Fz
            if i!=j: Fisher[j,i] = Fz
    return stats.FisherMatrix(Fisher,param_list)

def load_bao_experiment_rs_dV_diagonal(exp_name,data_path,boss_include=['6df','mgs','lowz','cmass']):
    if exp_name=='desi':
        zs,sig_pers = np.loadtxt(f'{data_path}/desi.txt',unpack=True)
    elif exp_name=='boss':
        array = np.genfromtxt(f'{data_path}/boss.txt',delimiter=',',dtype=None,encoding='utf-8')
        zs = []
        sig_pers = []
        for line in array:
            if line[0]  in boss_include: 
                zs.append(line[1])
                sig_pers.append(line[2])
    else:
        raise ValueError
    return np.asarray(zs),np.asarray(sig_pers)

def get_cls(params=None,lmax=3000,accurate=False,engine='camb',de='ppf',nonlinear=True):
    params = map_params(params,engine=engine)
    if engine=='camb':
        pars = set_camb_pars(params=params,de=de)
        if accurate:
            pars.set_accuracy(AccuracyBoost=3.0, lSampleBoost=1.0, lAccuracyBoost=3.0)
            pars.set_for_lmax(lmax=int(lmax+500), lens_potential_accuracy=3.0, max_eta_k=20000)
        else:
            pars.set_accuracy(AccuracyBoost=1.0, lSampleBoost=1.0, lAccuracyBoost=1.0)
            pars.set_for_lmax(lmax=int(lmax+500), lens_potential_accuracy=1, max_eta_k=2*lmax)
        if nonlinear:
            pars.NonLinear = model.NonLinear_both
        else:
            pars.NonLinear = model.NonLinear_none
        return load_theory(pars,lpad=lmax+2)
        
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
    """
    The matrix is input/output

    So it has elements like

    dAs/ds8, d(omc*h**2)/dom, d(omc*h**2)/ds8
    """

    if iparam==oparam: 
        return 1

    fs = fiducials
    h = fs['H0']/100.
    oc = fs['omch2'] / h**2.
    ob = fs['ombh2'] / h**2.
    if iparam=='omch2' and oparam=='om':
        return h**2.
    elif iparam=='ombh2' and oparam=='om':
        return h**2.
    elif iparam=='omch2' and oparam=='H0':
        return 2. * h * oc / 100.
    elif iparam=='ombh2' and oparam=='H0':
        return 2. * h * ob / 100.
    elif oparam=='s8':
        if iparam in ['As','ns','omch2','ombh2','ok','mnu','w0','wa','nnu','ctheta','thetastar','cs2']:
            val = np.loadtxt(derivpath+f"_fitderiv_s8_wrt_{iparam}.txt")
            assert val.size==1
            return 1./val.ravel()[0]
        elif iparam in ['tau','r']:
            return 0.
        else:
            raise ValueError
    

def reparameterize(Fmat,iparams,oparams,fiducials,deriv_path):
    """
    Re-parameterize a Fisher matrix.
    iparams must only contain CAMB primitives
    like As,ns,omch2,ombh2,H0,tau,mnu,w0,wa,omk,r.
    oparams can contain those as well as
    s8
    om
    """

    onum = len(oparams)
    inum = len(iparams)
    
    M = np.zeros((inum,onum))

    for i in range(inum):
        for j in range(onum):
            M[i,j] = get_trans_deriv(iparams[i],oparams[j],fiducials,deriv_path)
    print(M)
    return np.dot(np.dot(M,Fmat),M.T)


#def s8(params):
    

def load_theory(pars,lpad=9000):
    '''
    All ell and 2pi factors are stripped off.
    '''

    uSuffix = "unlensed_total"
    lSuffix = "total"

    results = camb.get_results(pars)
    cmbmat = results.get_cmb_power_spectra(pars,lmax=lpad,spectra=['total','unlensed_total','lens_potential'],raw_cl=True,CMB_unit='muK')

    theory = cosmology.TheorySpectra()
    for i,pol in enumerate(['TT','EE','BB','TE']):
        cls =cmbmat[lSuffix][:,i]
        ells = np.arange(0,len(cls),1)
        theory.loadCls(ells,cls,pol,lensed=True,interporder="linear",lpad=lpad,fill_zero=True)
        cls = cmbmat[uSuffix][:,i]
        ells = np.arange(0,len(cls),1)
        theory.loadCls(ells,cls,pol,lensed=False,interporder="linear",lpad=lpad,fill_zero=True)

    lensArr = cmbmat['lens_potential']
    cldd = lensArr[:,0]
    ells = np.arange(0,len(cldd),1)
    clkk = cldd.copy()
    clkk[1:] = clkk[1:] / 4. * (ells[1:]*(ells[1:]+1))**2.
    theory.loadGenericCls(ells,clkk,"kk",lpad=lpad,fill_zero=True)
    theory.dimensionless = False
    return theory





