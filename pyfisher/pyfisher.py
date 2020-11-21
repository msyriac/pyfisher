from __future__ import print_function
from orphics import cosmology,stats,maps,io
import numpy as np
import os,sys
import warnings
from pandas import DataFrame
import pandas as pd
import datetime


def prepare_output(args, message=""):
    output_path = args.output
    assert output_path.strip()[-1]!='/'
    io.mkdir(f'{output_path}')
    rname = os.path.basename(f'{output_path}')
    with open(f'{output_path}/info.log','w') as f:
        f.write(f'{message}\n')
        now = datetime.datetime.now()
        f.write(f'Current date and time : {now.strftime("%Y-%m-%d %H:%M:%S")}\n')
        for arg in vars(args):
            f.write(f'{arg} :  {getattr(args, arg)}\n')
        info = get_info(path=os.path.realpath(__file__))
        f.write(pretty_info(info))
    output_root = f'{output_path}/{rname}'
    return output_root

def get_saved_fisher(name,fsky=None,root_name='v20201120'):
    if name=='planck_lowell':
        fsky = 1 if fsky is None else fsky
        return fsky * read_fisher(f'{os.path.realpath(__file__)}/data/{root_name}_planck_low_ell_TT_fullsky.txt',delim=',')
    elif name=='planck_highell':
        fsky = 1 if fsky is None else fsky
        return fsky * read_fisher(f'{os.path.realpath(__file__)}/data/{root_name}_planck_high_ell_TTEETE_fullsky.txt',delim=',')
    elif name=='desi_bao':
        assert fsky is None
        return read_fisher(f'{os.path.realpath(__file__)}/data/{root_name}_desi_bao_fisher.txt',delim=',')
    elif name=='boss_bao':
        assert fsky is None
        return read_fisher(f'{os.path.realpath(__file__)}/data/{root_name}_boss_bao_fisher.txt',delim=',')




def check_fisher_sanity(fmat,param_list):
    Ny,Nx = fmat.shape
    assert Ny==Nx
    assert Ny==len(param_list)
    assert len(param_list)==len(set(param_list))

def write_fisher(filename,fmat,delim=','):
    np.savetxt(filename,fmat,header=(delim).join(fmat.params),delimiter=delim)

def read_fisher(csv_file,delimiter=','):
    fmat = np.loadtxt(csv_file,delimiter=delimiter)
    with open(csv_file) as f:
        fline = f.readline()
    fline = fline.replace("#","")
    columns = fline.strip().split(delimiter)
    assert len(set(columns)) == len(columns)
    return FisherMatrix(fmat = fmat,param_list = columns)

def rename_fisher(fmat,pmapping):
    old_params = fmat.params
    new_params = list(old_params)
    for key in pmapping.keys():
        if key not in old_params: continue
        i = old_params.index(key)
        new_params[i] = pmapping[key]
    return FisherMatrix(fmat=fmat.values,param_list=new_params)
    
class FisherMatrix(DataFrame):
    """
    A Fisher Matrix object that subclasses pandas.DataFrame.
    This is essentially just a structured array that
    has identical column and row labels.

    You can initialize an empty one like:
    >> params = ['H0','om','sigma8']
    >> F = FisherMatrix(np.zeros((len(params),len(params))),params)
    
    where params is a list of parameter names. If you already have a
    Fisher matrix 'Fmatrix' whose diagonal parameter order is specified by
    the list 'params', then you can initialize this object as:
    
    >> F = FisherMatrix(Fmatrix,params)
    
    This makes the code 'aware' of the parameter order in a way that makes
    handling combinations of Fishers a lot easier.
    
    You can set individual elements like:
    
    >> F['s8']['H0'] = 1.

    Once you've populated the entries, you can do things like:
    >> Ftot = F1 + F2
    i.e. add Fisher matrices. The nice property here is that you needn't
    have initialized these objects with the same list of parameters!
    They can, for example, have mutually exclusive parameters, in which
    case you end up with some reordering of a block diagonal Fisher matrix.
    In the general case, of two overlapping parameter lists that don't
    have the same ordering, pandas will make sure the objects are added
    correctly.

    WARNING: No other operation other than addition and multiplication is overloaded. Subtraction
    for instance will give unpredictable behaviour. (Will likely introduce
    NaNs) But you shouldn't be subtracting Fisher matrices anyway!

    You can add a gaussian prior to a parameter:
    >> F.add_prior('H0',2.0)

    You can drop an entire parameter (which removes that row and column):
    >> F.delete('s8')
    which does it in place.

    If you want to preserve the original before modifying, you can
    >> Forig = F.copy()

    You can get marginalized errors on each parameter as a dict:
    >> sigmas = F.sigmas()


    """
    
    def __init__(self,fmat,param_list,delete_params=None,prior_dict=None):
        """
        fmat            -- (n,n) shape numpy array containing initial Fisher matrix for n parameters
        param_list      -- n-element list specifying diagonal order of fmat
        delete_params   -- list of names of parameters you would like to delete from this 
                        Fisher matrix when you initialize it. 
        prior_dict      -- a dictionary that maps names of parameters to 1-sigma prior values
                        you would like to add on initialization. This can also be done later with the 
                        add_prior function.
	"""
	
	
        check_fisher_sanity(fmat,param_list)
        pd.DataFrame.__init__(self,fmat.copy(),columns=param_list,index=param_list)
        try:
            a = self.params
            raise ValueError # self.params should not already exist
        except:
            pass
        self.params = param_list
            
        cols = self.columns.tolist()
        ind = self.index.tolist()
        assert set(self.params)==set(cols)
        assert set(self.params)==set(ind)

        if delete_params is not None:
            self.delete(delete_params)
        if prior_dict is not None:
            for prior in prior_dict.keys():
                self.add_prior(prior,prior_dict[prior])

            
    def copy(self):
        """
        >> Fnew = F.copy()
        will create an independent Fnew that is not a view of the original.
        """
        f = FisherMatrix(pd.DataFrame.copy(self), list(self.params))
        return f

    def __radd__(self,other):
        return self._add(other,radd=True)

    def __add__(self,other):
        return self._add(other,radd=False)

    def __mul__(self,other):
        return FisherMatrix(self.values*other,self.columns.tolist())

    def __rmul__(self,other):
        return FisherMatrix(self.values*other,self.columns.tolist())

    def _add(self,other,radd=False):
        if other is None: return self
        F1 = pd.DataFrame(self.values,columns=self.params,index=self.params)
        F2 = pd.DataFrame(other.values,columns=other.params,index=other.params)
        if radd:
            new_fpd = pd.DataFrame.radd(F1,F2,fill_value=0)
        else:
            new_fpd = pd.DataFrame.add(F1,F2,fill_value=0)
        return FisherMatrix(new_fpd.values,new_fpd.columns.tolist())

    def add_prior(self,param,prior):
        """
        Adds 1-sigma value 'prior' to the parameter name specified by 'param'
        """
        self[param][param] += 1./prior**2.
        
    def sigmas(self):
        """
        Returns marginalized 1-sigma uncertainties on each parameter in the Fisher matrix.
        """
        finv = np.linalg.inv(self.values)
        errs = np.diagonal(finv)**(0.5)
        return dict(zip(self.params,errs))
    
    def delete(self,params):
        """
        Given a list of parameter names 'params', deletes these from the Fisher matrix.
        """
        self.drop(labels=params,axis=0,inplace=True)
        self.drop(labels=params,axis=1,inplace=True)
        self.params = self.columns.tolist()
        assert set(self.index.tolist())==set(self.params)

    def reordered(self,params):
        # Return a reordered version of self
        return self[params].T[params]

    def marge_var_2param(self,param1,param2):
        """
        Returns the sub-matrix corresponding to two parameters param1 and param2.
        Useful for contour plots.
        """
        finv = np.linalg.inv(self.values)
        i = self.params.index(param1)
        j = self.params.index(param2)
        chi211 = finv[i,i]
        chi222 = finv[j,j]
        chi212 = finv[i,j]
        return np.array([[chi211,chi212],[chi212,chi222]])


def get_planck_cmb_fisher(param_list,bin_edges,specs,root_name,fsky,interpolate=True):
    ells = np.arange(0,bin_edges.max()+1)
    nls = get_planck_nls(ells)
    cls = load_theory_dict(f'{root_name}_fiducial.txt',ells)
    dcls = load_derivs(root_name,param_list,ells)
    return band_fisher(param_list,bin_edges,specs,cls,nls,dcls,interpolate=interpolate)  * fsky


def load_derivs(root_name,param_list,ells):
    dcls = {}
    for param in param_list:
        dcls[param] = load_theory_dict(f'{root_name}_{param}_deriv.txt',ells)
    return dcls


def load_theory_dict(fname,ells):
    cls = {}
    ells,tt,ee,bb,te,kk = np.loadtxt(fname,unpack=True)
    cls['TT'] = maps.interp(ells,tt)
    cls['EE'] = maps.interp(ells,ee)
    cls['BB'] = maps.interp(ells,bb)
    cls['TE'] = maps.interp(ells,te)
    cls['kk'] = maps.interp(ells,kk)
    return cls


def get_planck_nls(ells):
    beams_T =  [33.,23.,14.,10.,7.,5.,5.]
    uK_arcmins_T = [145.,149.,137.,65.,43.,66.,200.]
    beams_P =  [14.,10.,7.,5.,5.]
    uK_arcmins_P = [450.,103.,81.,134.,406.]
    Ns_TT = np.asarray([(uK_arcmin*np.pi/180./60.)**2./maps.gauss_beam(ells,fwhm)**2. for uK_arcmin,fwhm in zip(uK_arcmins_T,beams_T)])
    Ns_PP = np.asarray([(uK_arcmin*np.pi/180./60.)**2./maps.gauss_beam(ells,fwhm)**2. for uK_arcmin,fwhm in zip(uK_arcmins_P,beams_P)])
    N_TT = 1./(1./Ns_TT).sum(axis=0)
    N_PP = 1./(1./Ns_PP).sum(axis=0)
    nls = {}
    nls['TT'] = maps.interp(ells,N_TT)
    nls['EE'] = maps.interp(ells,N_PP)
    nls['BB'] = maps.interp(ells,N_PP)
    return nls


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

    return FisherMatrix(Fisher,param_list)




def get_param_info(param_file,exclude=None):
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
    import camb
    from camb import model
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
    import camb
    from camb import model
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
        from classy import Class
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
    return FisherMatrix(Fisher,param_list)

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
        from classy import Class
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

    import camb
    from camb import model
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


# COPIED FROM ACTSIMS
def pretty_info(info):
    name = info['package'] if info['package'] is not None else info['path']
    pstr = f'\n{name}'
    pstr = pstr + '\n'+''.join(["=" for x in range(len(name))])
    for key in info.keys():
        if key=='package': continue
        pstr = pstr + '\n' + f'\t{key:<10}{str(info[key]):<40}'
    return pstr

# COPIED FROM ACTSIMS
def get_info(package=None,path=None,validate=True):
    import git
    import importlib
    info = {}
    if package is None:
        assert path is not None, "One of package or path must be specified."
        path = os.path.dirname(path)
        version = None
    else:
        mod = importlib.import_module(package)
        try:
            version = mod.__version__
        except AttributeError:
            version = None
        path = mod.__file__
        path = os.path.dirname(path)
    info['package'] = package
    info['path'] = path
    info['version'] = version
    try:
        repo = git.Repo(path,search_parent_directories=True)
        is_git = True
    except git.exc.InvalidGitRepositoryError:
        is_git = False
    info['is_git'] = is_git
    if is_git:
        chash = str(repo.head.commit)
        untracked = len(repo.untracked_files)>0
        changes = len(repo.index.diff(None))>0
        branch = str(repo.active_branch)
        info['hash'] = chash
        info['untracked'] = untracked
        info['changes'] = changes
        info['branch'] = branch
    else:
        if validate:
            assert version is not None
            assert 'site-packages' in path
    return info
    