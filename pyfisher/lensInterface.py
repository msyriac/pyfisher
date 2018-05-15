import numpy as np
from scipy.interpolate import interp1d
import sys
from orphics.io import list_from_config

def lensNoise(Config,expName,lensName,beamOverride=None,noiseTOverride=None,lkneeTOverride=None,lkneePOverride=None,alphaTOverride=None,alphaPOverride=None,tellminOverride=None,pellminOverride=None,tellmaxOverride=None,pellmaxOverride=None,deg=5.,px=1.0,gradCut=10000,bigell=9000,plot=False,theoryOverride=None,lensedEqualsUnlensed=True,noiseFuncT=None,noiseFuncP=None,nIter=np.inf):

    from orphics.io import list_from_config

    beam = list_from_config(Config,expName,'beams')
    noise = list_from_config(Config,expName,'noises')
    freq = list_from_config(Config,expName,'freqs')
    lkneeT,lkneeP = list_from_config(Config,expName,'lknee')
    alphaT,alphaP = list_from_config(Config,expName,'alpha')
    if (noiseFuncT is None) and (noiseFuncP is None):
        print 'Not using noise files for generating lensing noise'
    else:
        print 'Using noise files for generating lensing noise'
    tellmin,tellmax = list_from_config(Config,expName,'tellrange')
    if tellminOverride is not None: tellmin = tellminOverride
    if tellmaxOverride is not None: tellmax = tellmaxOverride
    pellmin,pellmax = list_from_config(Config,expName,'pellrange')
    if pellminOverride is not None: pellmin = pellminOverride
    if pellmaxOverride is not None: pellmax = pellmaxOverride
    lmax = int(Config.getfloat(expName,'lmax'))

    pols = Config.get(lensName,'polList').split(',')
    freq_to_use = Config.getfloat(lensName,'freq')
    ind = np.where(np.isclose(freq,freq_to_use))
    beamFind = np.array(beam)[ind]
    noiseFind = np.array(noise)[ind]
    assert beamFind.size==1
    assert noiseFind.size==1
    if beamOverride is not None:
        beamX = beamY = beamOverride
    else:
        beamX = beamY = beamFind[0]
    if noiseTOverride is not None:
        noiseTX = noiseTY = noiseTOverride
    else:
        noiseTX = noiseTY = noiseFind[0]
    if lkneeTOverride is not None: lkneeT = lkneeTOverride
    if lkneePOverride is not None: lkneeP = lkneePOverride
    if alphaTOverride is not None: alphaT = alphaTOverride
    if alphaPOverride is not None: alphaP = alphaPOverride

    from orphics.lensing import NlGenerator,getMax
    from orphics import maps

    deg = 20
    dell = 10
    kellmin = 2
    kellmax = max(tellmax,pellmax)

    shape,wcs = maps.rect_geometry(width_deg = deg, px_res_arcmin=px)
    eff_at = 60.
    
    kmin,kmax = list_from_config(Config,lensName,'Lrange')

    #print kellmin,kellmax,pellmin,pellmax,kmin,kmax
    
    if theoryOverride is None:
        from orphics.cosmology import Cosmology
        cc = Cosmology(lmax=int(kellmax),pickling=True)
        theory = cc.theory
    else:
        theory = theoryOverride
        cc = None
    bin_edges = np.arange(kellmin,kellmax,dell)
    myNls = NlGenerator(shape,wcs,theory,bin_edges,gradCut=gradCut,bigell=bigell,unlensedEqualsLensed=lensedEqualsUnlensed)
    myNls.updateNoise(beamX,noiseTX,np.sqrt(2.)*noiseTX,tellmin,tellmax,pellmin,pellmax,beamY=beamY,noiseTY=noiseTY,noisePY=np.sqrt(2.)*noiseTY,lkneesX=(lkneeT,lkneeP),lkneesY=(lkneeT,lkneeP),alphasX=(alphaT,alphaP),alphasY=(alphaT,alphaP),noiseFuncTX=noiseFuncT,noiseFuncTY=noiseFuncT,noiseFuncPX=noiseFuncP,noiseFuncPY=noiseFuncP)

    
    lsmv,Nlmv,ells,dclbb,efficiency = myNls.getNlIterative(pols,pellmin,pellmax,dell=dell,halo=True,plot=plot,max_iterations=nIter,eff_at=eff_at,kappa_min=kmin,kappa_max=kmax)

     
    return lsmv,Nlmv,ells,dclbb,efficiency,cc


