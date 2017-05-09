import os
import itertools
import numpy as np
import sys
from scipy.interpolate import interp1d
import ConfigParser
import traceback
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
from scipy.linalg import block_diag

def rSigma(fsky,ellBBRange,fnBB,dCls,lclbb,ClBB=lambda x: 0.):

    # lclbb, ell=2 starts at 0
    # dCls, ClBB, ell=2 starts at 2
    

    Fisher = 0.
    for ell in ellBBRange:

        dBBdr = dCls(ell)
        deltaCl = ClBB(ell) + fnBB(ell) + lclbb(ell)

                
        if dBBdr==0.:
            Frrell = 0.
        else:    
            Frrell = fsky*(2.*ell+1.)*(dBBdr**2.)/2./(deltaCl**2.)
    
        Fisher += Frrell

    return np.sqrt(1./Fisher)


def CovFromVecs(Cls,ell,nTT=0.,nEE=0.,nkk=0.,lensing=False):
    '''
    For calculating Fisher through Eq.4 of 1402.4108 ("compact" option)
    Pros: easily extendable to optical cross-correlations
    Cons: doesn't currently reproduce non-compact option exactly
    '''

    TT = Cls[ell,0] + nTT
    EE = Cls[ell,1] + nEE
    TE = Cls[ell,3]


    kk = (Cls[ell,4] + nkk)
    kt = (Cls[ell,5])

    if lensing:
        mat = np.array([[TT,TE,kt],
                       [TE,EE,0.],
                       [kt,0.,kk]])
    else:
        mat = np.array([[TT,TE],
                       [TE,EE]])
    return mat



def calcFisher(paramList,ellrange,fidCls,dCls,fnTT,fnEE,fnKK,fsky,lensing=True,verbose=True):
    numParams = len(paramList)
    
    Cls = []
    nCls = []
    # Loop through each unique parameter combination
    Fisher = np.zeros((numParams,numParams))
    paramCombs = itertools.combinations_with_replacement(paramList,2)
    for param1,param2 in paramCombs:
        i = paramList.index(param1)
        j = paramList.index(param2)
        Fell = 0.
        for ell in ellrange:

            nTT = fnTT(ell)
            nEE = fnEE(ell)

            nkk = fnKK(ell)
            Cov = CovFromVecs(fidCls,ell,nTT,nEE,nkk=nkk,lensing=lensing)
            dCov1 = CovFromVecs(dCls[param1],ell,lensing=lensing)
            dCov2 = CovFromVecs(dCls[param2],ell,lensing=lensing)
            InvCov = np.linalg.inv(Cov)
            Fell += (2.*ell+1.) * fsky * np.trace(np.dot(np.dot(InvCov,dCov1),np.dot(InvCov,dCov2))) /2.


        Fisher[i,j] = Fell
        Fisher[j,i] = Fell


    return Fisher


def calcBvec(paramList,ellrange,sigCls,fidCls,dCls,fnTT,fnEE,fnKK,fsky,lensing=True,verbose=True):
    numParams = len(paramList)
    
    Cls = []
    nCls = []
    # Loop through each unique parameter combination
    Bvec = np.zeros((numParams,1))
    for param in paramList:
        i = paramList.index(param )
        B= 0.
        for ell in ellrange:

            nTT = fnTT(ell)
            nEE = fnEE(ell)

            nkk = fnKK(ell)
            Cov = CovFromVecs(fidCls,ell,nTT,nEE,nkk=nkk,lensing=lensing)
            dCov1 = CovFromVecs(sigCls,ell,lensing=lensing)
            dCov2 = CovFromVecs(dCls[param],ell,lensing=lensing)
            InvCov = np.linalg.inv(Cov)
            B += (2.*ell+1.) * fsky * np.trace(np.dot(np.dot(InvCov,dCov1),np.dot(InvCov,dCov2))) /2.


        Bvec[i] = B
    return Bvec

def calcBias(Fisher,Bvec):
    Finv = np.linalg.inv(Fisher)
    return np.dot(Finv,Bvec)


def tryLoad(filepath,delimiter=None):
    try:
        return np.loadtxt(filepath,delimiter=delimiter)
    except:
        try:
            return np.loadtxt('output/'+filepath,delimiter=delimiter)
        except:
            return np.loadtxt(os.environ['FISHER_DIR']+'/output/'+filepath,delimiter=delimiter)
def loadFishers(filepaths):
    totFisher = 0.
    for filepath in filepaths:
        F = np.loadtxt(filepath)
        totFisher += F

    return totFisher

def noiseFromConfig(Config,expName,TCMB=2.7255e6,beamsOverride=None,noisesOverride=None,lkneeTOverride=None,lkneePOverride=None,alphaTOverride=None,alphaPOverride=None,tellminOverride=None,pellminOverride=None,tellmaxOverride=None,pellmaxOverride=None):
    import orphics.tools.cmb as cmb
    from orphics.tools.io import dictFromSection, listFromConfig

    tellmin,tellmax = listFromConfig(Config,expName,'tellrange')
    if tellminOverride is not None: tellmin = tellminOverride
    if tellmaxOverride is not None: tellmax = tellmaxOverride
    pellmin,pellmax = listFromConfig(Config,expName,'pellrange')
    if pellminOverride is not None: pellmin = pellminOverride
    if pellmaxOverride is not None: pellmax = pellmaxOverride
    if beamsOverride is not None:
        beams = beamsOverride
    else:
        beams = listFromConfig(Config,expName,'beams')
    if noisesOverride is not None:
        noises = noisesOverride
    else:
        noises = listFromConfig(Config,expName,'noises')
    lkneeT,lkneeP = listFromConfig(Config,expName,'lknee')
    alphaT,alphaP = listFromConfig(Config,expName,'alpha')
    if lkneeTOverride is not None: lkneeT = lkneeTOverride
    if lkneePOverride is not None: lkneeP = lkneePOverride
    if alphaTOverride is not None: alphaT = alphaTOverride
    if alphaPOverride is not None: alphaP = alphaPOverride

    invnTTs = 0.
    invnEEs = 0.
    for beam,noise in zip(beams,noises):
       invnTTs += 1./cmb.noise_func(np.arange(tellmin,tellmax),beam,noise,lknee=lkneeT,alpha=alphaT,TCMB=TCMB)
       invnEEs += 1./cmb.noise_func(np.arange(pellmin,pellmax),beam,noise*np.sqrt(2.),lknee=lkneeP,alpha=alphaP,TCMB=TCMB)

    fnTT = interp1d(np.arange(tellmin,tellmax),1./invnTTs,bounds_error=False,fill_value=np.inf)
    fnEE = interp1d(np.arange(pellmin,pellmax),1./invnEEs,bounds_error=False,fill_value=np.inf)
    
    return fnTT, fnEE

def testAgainstKSmith(pmax,beamFWHMArcmin,dCls,lclbb,rExp,rInFid,fCls,fsky):
    from orphics.tools.io import Plotter
    pl = Plotter(scaleX='log',scaleY='log')
    pnoiserange = np.logspace(np.log10(0.5),np.log10(50.),num=100)
    for pmin in [2,5,10,40]:
        sigs = []
        for deltaP in pnoiserange:
            ellBBRange = range(pmin,pmax)

            sigs.append(rSigma(fsky,ellBBRange,beamFWHMArcmin,deltaP,dCls[:,2],lclbb,rExp*fCls[:,2]/rInFid))

        kn, kr = np.loadtxt("data/k"+str(pmin)+".csv",delimiter=',',unpack=True)
        pl.add(kn,kr,ls='--')
        pl.add(pnoiserange,sigs,label="$\\ell_{\mathrm{min}}="+str(pmin)+"$")

    pl.legendOn()
    pl._ax.set_xlim(0.5,50.)
    pl._ax.set_ylim(1.e-5,1.e-1)
    pl.done("kenplot.png")    
