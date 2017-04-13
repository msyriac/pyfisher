import numpy as np
from scipy.interpolate import interp1d

def lensNoise(Config,expName,lensName,beamOverride=None,noiseTOverride=None,lkneeTOverride=None,lkneePOverride=None,alphaTOverride=None,alphaPOverride=None):

    from orphics.tools.io import dictFromSection, listFromConfig

    beam = listFromConfig(Config,expName,'beams')
    noise = listFromConfig(Config,expName,'noises')
    freq = listFromConfig(Config,expName,'freqs')
    lkneeT,lkneeP = listFromConfig(Config,expName,'lknee')
    alphaT,alphaP = listFromConfig(Config,expName,'alpha')
    tellmin,tellmax = listFromConfig(Config,expName,'tellrange')
    pellmin,pellmax = listFromConfig(Config,expName,'pellrange')
    lmax = int(Config.getfloat(expName,'lmax'))

    pols = Config.get(lensName,'polList').split(',')
    delens = Config.getboolean(lensName,'delens')
    freq_to_use = Config.getfloat(lensName,'freq')
    ind = np.where(np.isclose(freq,freq_to_use))
    beamFind = np.array(beam)[ind]
    noiseFind = np.array(noise)[ind]
    assert beamFind.size==1
    assert noiseFind.size==1
    if beamOverride is not None:
        beamX = beamY = beam
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

    from orphics.tools.cmb import loadTheorySpectraFromCAMB
    import flipper.liteMap as lm
    from alhazen.quadraticEstimator import NlGenerator,getMax
    deg = 10.
    px = 0.5
    dell = 20
    gradCut = 10000
    kmin = 100
    # cambRoot = "data/ell28k_highacc"
    # theory = loadTheorySpectraFromCAMB(cambRoot,unlensedEqualsLensed=False,useTotal=False,lpad=9000)
    from orphics.theory.cosmology import Cosmology
    cc = Cosmology(lmax=6000,pickling=True)
    theory = cc.theory
    lmap = lm.makeEmptyCEATemplate(raSizeDeg=deg, decSizeDeg=deg,pixScaleXarcmin=px,pixScaleYarcmin=px)

    Nleach = {}
    kmaxes = []
    for polComb in pols:
        kmax = getMax(polComb,tellmax,pellmax)
        bin_edges = np.arange(kmin,kmax,dell)+dell

        myNls = NlGenerator(lmap,theory,bin_edges,gradCut=gradCut)
        myNls.updateNoise(beamX,noiseTX,np.sqrt(2.)*noiseTX,tellmin,tellmax,pellmin,pellmax,beamY=beamY,noiseTY=noiseTY,noisePY=np.sqrt(2.)*noiseTY,lkneesX=(lkneeT,lkneeP),lkneesY=(lkneeT,lkneeP),alphasX=(alphaT,alphaP),alphasY=(alphaT,alphaP))

        if (polComb=='EB' or polComb=='TB') and (delens):
            ls, Nls, eff = myNls.iterativeDelens(polComb,1.0,True)
        else:
            ls,Nls = myNls.getNl(polComb=polComb,halo=True)

        Nleach[polComb] = (ls,Nls)
        kmaxes.append(kmax)

    bin_edges = np.arange(kmin,max(kmaxes),dell)+dell
    Nlmvinv = 0.
    from scipy.interpolate import interp1d

    # from orphics.tools.io import Plotter
    # ellkk = np.arange(2,9000,1)
    # Clkk = theory.gCl("kk",ellkk)    
    # pl = Plotter(scaleY='log',scaleX='log')
    # pl.add(ellkk,4.*Clkk/2./np.pi)


    for polComb in pols:
        ls,Nls = Nleach[polComb]
        nlfunc = interp1d(ls,Nls,bounds_error=False,fill_value=np.inf)
        Nleval = nlfunc(bin_edges)
        Nlmvinv += np.nan_to_num(1./Nleval)
        # pl.add(ls,4.*Nls/2./np.pi,label=polComb)

    Nlmv = np.nan_to_num(1./Nlmvinv)
    ls = bin_edges[1:-1]
    Nls = Nlmv[1:-1]

    # from orphics.tools.stats import bin1D
    # binner1d = bin1D(bin_edges)
    # ellcls , clkk_binned = binner1d.binned(ellkk,Clkk)



    # pl.add(ellcls,4.*clkk_binned/2./np.pi,ls="none",marker="x")
    # pl.add(ellcls,4.*clkk_binned/2./np.pi,ls="none",marker="x")
    # pl.add(ls,4.*Nls/2./np.pi,ls="--")
    # np.savetxt("output/nlsave_"+expName+"_"+lensName+".txt",np.vstack((ls,Nls)).transpose())

    # ####Nls += clkk_binned[:-1]
    # np.savetxt("output/nlsaveTot_"+expName+"_"+lensName+".txt",np.vstack((ls,Nls)).transpose())
    # pl.add(ls,4.*Nls/2./np.pi,ls="-.")

    # pl.legendOn(loc='lower left',labsize=10)
    # pl.done("output/Nl_"+expName+lensName+".png")
    
    return ls,Nls
