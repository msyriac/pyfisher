'''
Output BAO Fisher matrix into .csv file in output/
'''
import camb
from camb import model, initialpower
import numpy as np
import sys
import ConfigParser
import itertools

from pyfisher.core import getHubbleCosmology

def getBAOCamb(zrange,params,AccuracyBoost=False):
    '''
    Get BAO's fk=rs/D_V from CAMB. Return 1D array of fk for all redshifts z
    '''
    pars = camb.CAMBparams()
    if AccuracyBoost:
        pars.set_accuracy(AccuracyBoost=2.0, lSampleBoost=2.0, lAccuracyBoost=2.0)
        #pars.set_accuracy(AccuracyBoost=3.0, lSampleBoost=3.0, lAccuracyBoost=3.0)

    pars.set_dark_energy(w=params['w'],wa=params['wa'],dark_energy_model='ppf')
    try:
        theta = params['theta100']/100.
        H0 = None
    except:
        H0 = params['H0']
        theta = None
    pars.set_cosmology(H0=H0, cosmomc_theta=theta, ombh2=params['ombh2'], omch2=params['omch2'], tau=params['tau'],mnu=params['mnu'],nnu=params['nnu'],omk=params['omk'],num_massive_neutrinos=int(params['num_massive_neutrinos']),TCMB=params['TCMB'])
    pars.InitPower.set_params(As=params['As'],ns=params['ns'],r=params['r'])
    #pars.set_for_lmax(lmax=int(params['lmax']), lens_potential_accuracy=1, max_eta_k=2*params['lmax'])
    camb.set_z_outputs(zrange)
    results = camb.get_results(pars)
    fk = results.get_background_outputs() #results.get_BAO(zrange,pars)[:,0]
    #print fk
    fk = np.nan_to_num(fk[:,0])
    return fk

def main(argv):
    verbose=True

    # Read Config
    iniFile = "input/" + argv[0] #"makeDefaultsBAO_szar.ini"
    Config = ConfigParser.SafeConfigParser()
    Config.optionxform = str
    Config.read(iniFile)
    out_pre = Config.get('general','output_prefix')
    AccuracyBoost = Config.getboolean('general','AccuracyBoost')
    paramList = []
    fparams = {}
    stepSizes = {}
    for (key, val) in Config.items('camb'):
        if ',' in val:
            param, step = val.split(',')
            paramList.append(key)
            fparams[key] = float(param)
            stepSizes[key] = float(step)
        else:
            fparams[key] = float(val)    
    zrange = [float(x) for x in Config.get('DESI','redshift').split(',')]
    sigmafk = [float(x)/1000. for x in Config.get('DESI','sigmafkx1000').split(',')]
	
    # Uncomment if need to calculate Derivs

    # Save fiducials
    print "Calculating and saving fiducial cosmology..."
    # if not('H0' in fparams):
    #     fparams['H0'] = getHubbleCosmology(theta=fparams['theta100'],params=fparams)
    fidfk = getBAOCamb(zrange,fparams,AccuracyBoost=AccuracyBoost)
    np.savetxt("output/"+out_pre+"_fidfk.csv",fidfk,delimiter=",")


    # Calculate and save derivatives
    for paramName in paramList:
        h = stepSizes[paramName]        
        print "Calculating forward difference for ", paramName
        pparams = fparams.copy()
        pparams[paramName] = fparams[paramName] + 0.5*h
        # if paramName=='theta100':
        #     pparams['H0'] = getHubbleCosmology(theta=pparams['theta100'],params=pparams)
        pfk = getBAOCamb(zrange,pparams,AccuracyBoost=AccuracyBoost)
    
    
        print "Calculating backward difference for ", paramName
        mparams = fparams.copy()
        mparams[paramName] = fparams[paramName] - 0.5*h
        # if paramName=='theta100':
        #     mparams['H0'] = getHubbleCosmology(theta=mparams['theta100'],params=mparams)
        mfk = getBAOCamb(zrange,mparams,AccuracyBoost=AccuracyBoost)
    
        dfk = (pfk-mfk)/h
    
        np.savetxt("output/"+out_pre+"_dfk_"+paramName+".csv",dfk,delimiter=",")
        	
    # Calculate Fisher Matrix
    paramCombs = itertools.combinations_with_replacement(paramList,2)
    Fisher = np.zeros((len(paramList),len(paramList)))
    dfks = {}
    for paramName in paramList:
        dfks[paramName] = np.loadtxt("output/"+out_pre+"_dfk_"+paramName+".csv",delimiter=",")
        
    for param1,param2 in paramCombs:
        if verbose: print "Parameter combination : ", param1,param2
        i = paramList.index(param1)
        j = paramList.index(param2)
        Fz = 0.		
        for k in range(0,len(zrange)):
            dfk1 = dfks[param1][k]
            dfk2 = dfks[param2][k]
            Fz += dfk1*dfk2/sigmafk[k]**2.
            #if verbose: print "dfk1,dfk2,sigmafk,Fz:",dfk1*dfk2/sigmafk[k]**2,Fz
        Fisher[i,j] = Fz
        Fisher[j,i] = Fz
    if verbose:
        print Fisher
        print np.diagonal(Fisher)
    np.savetxt("output/"+out_pre+"_Fisher.csv",Fisher,delimiter=",")
    
if (__name__ == "__main__"):
    main(sys.argv[1:])
