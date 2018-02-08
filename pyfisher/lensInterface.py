import numpy as np
from scipy.interpolate import interp1d


def lensNoise(Config,expName,lensName,beamOverride=None,noiseTOverride=None,lkneeTOverride=None,lkneePOverride=None,alphaTOverride=None,alphaPOverride=None,tellminOverride=None,pellminOverride=None,tellmaxOverride=None,pellmaxOverride=None,px=1.0,gradCut=10000,bigell=9000,plot=False):

    from orphics.tools.io import dictFromSection, listFromConfig

    beam = listFromConfig(Config,expName,'beams')
    noise = listFromConfig(Config,expName,'noises')
    freq = listFromConfig(Config,expName,'freqs')
    lkneeT,lkneeP = listFromConfig(Config,expName,'lknee')
    alphaT,alphaP = listFromConfig(Config,expName,'alpha')
    tellmin,tellmax = listFromConfig(Config,expName,'tellrange')
    if tellminOverride is not None: tellmin = tellminOverride
    if tellmaxOverride is not None: tellmax = tellmaxOverride
    pellmin,pellmax = listFromConfig(Config,expName,'pellrange')
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

    import flipper.liteMap as lm
    from alhazen.quadraticEstimator import NlGenerator,getMax
    deg = 5.
    #px = 1.0
    dell = 10
    #gradCut = 10000
    kellmin = 10
    shape,wcs = maps.rect_geometry(width_deg = deg, px_res_arcmin=px)
    kellmax = max(tellmax,pellmax)
    from orphics.theory.cosmology import Cosmology
    cc = Cosmology(lmax=int(kellmax),pickling=True)
    theory = cc.theory
    bin_edges = np.arange(kellmin,kellmax,dell)
    myNls = NlGenerator(lmap,theory,bin_edges,gradCut=gradCut,bigell=bigell)
    myNls.updateNoise(beamX,noiseTX,np.sqrt(2.)*noiseTX,tellmin,tellmax,pellmin,pellmax,beamY=beamY,noiseTY=noiseTY,noisePY=np.sqrt(2.)*noiseTY,lkneesX=(lkneeT,lkneeP),lkneesY=(lkneeT,lkneeP),alphasX=(alphaT,alphaP),alphasY=(alphaT,alphaP))

    lsmv,Nlmv,ells,dclbb,efficiency = myNls.getNlIterative(pols,kellmin,kellmax,tellmax,pellmin,pellmax,dell=dell,halo=True,plot=plot)

     
    return lsmv,Nlmv,ells,dclbb,efficiency,cc


