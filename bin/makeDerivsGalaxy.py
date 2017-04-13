'''
Output galaxy Fisher matrix into .csv file in output/
'''
import numpy as np
import sys,os
import ConfigParser
import itertools

#from cambRedCall import cambRed
from cambCall import cambInterface

def getClsCambRed(out_pre,templateIni,params,AccuracyBoost=False,savefid=False,seed=0):
    option = 1
    cRed = cambInterface(out_pre,templateIni,cambRoot=os.environ['CAMBREDNN_DIR'],option=option,seed=seed)
    for key in params:
        cRed.setParam(key,params[key])
    #sys.exit()
    cRed.call(suppress=False)
    #print 'call cambRed done!'
    Cls = cRed.getCls()
    if savefid:
        os.system('cp $CAMBREDNN_DIR/'+out_pre+'_params.ini $FISHER_DIR/output/'+out_pre+'_galaxy_fid.csv')
    cRed.done()
    return Cls

def main(argv):
    verbose=True
    nameMap={'H0':'hubble','YHe':'helium_fraction','nnu':'massless_neutrinos','s_pivot':'pivot_scalar','t_pivot':'pivot_tensor','As':'scalar_amp(1)','ns':'scalar_spectral_index(1)','tau':'re_optical_depth','num_massive_neutrinos':'massive_neutrinos','TCMB':'temp_cmb','lmax':'l_max_scalar'}
    inv_nameMap = {v: k for k, v in nameMap.items()}

    # Read Config
    iniFile = "input/makeDerivsGalaxy.ini"
    Config = ConfigParser.SafeConfigParser()
    Config.optionxform = str
    Config.read(iniFile)
    out_pre = Config.get('general','output_prefix')
    AccuracyBoost = Config.getboolean('general','AccuracyBoost')
    templateIni = os.environ['CAMBREDNN_DIR']+Config.get('general','templateIni')
    paramList = []
    fparams = {}
    stepSizes = {}
    for (key, val) in Config.items('cambRed'):
        if key in nameMap:
            key = nameMap[key]
        if ',' in val:
            param, step = val.split(',')
            paramList.append(key)
            fparams[key] = float(param)
            stepSizes[key] = float(step)
        else:
            if key == 'l_max_scalar':
                fparams[key] = int(val)
            else:
                fparams[key] = float(val)    

    # Uncomment if need to calculate Derivs
    print fparams
    print paramList
    # Save fiducials
    print "Calculating and saving fiducial cosmology..."
    if not('omnuh2' in fparams) and ('mnu' in fparams):
        fparams['omnuh2'] = round(fparams['mnu']/93.14,6)
    fCls = getClsCambRed(out_pre,templateIni,fparams,AccuracyBoost=AccuracyBoost,savefid=True)
    np.savetxt("output/"+out_pre+"_galaxy_fCls.csv",fCls,delimiter=",")

    # Calculate and save derivatives
    for paramName in paramList:
        h = stepSizes[paramName]        
        print "Calculating forward difference for ", paramName
        pparams = fparams.copy()
        pparams[paramName] = fparams[paramName] + 0.5*h
        if paramName =='mnu':
            pparams['omnuh2'] = round(pparams['mnu']/93.14,6)
        pCls = getClsCambRed(out_pre,templateIni,pparams,AccuracyBoost=AccuracyBoost)
    
    
        print "Calculating backward difference for ", paramName
        mparams = fparams.copy()
        mparams[paramName] = fparams[paramName] - 0.5*h
        if paramName =='mnu':
            mparams['omnuh2'] = round(mparams['mnu']/93.14,6)
        mCls = getClsCambRed(out_pre,templateIni,mparams,AccuracyBoost=AccuracyBoost)

        dCls = (pCls-mCls)/h
        if paramName in inv_nameMap:
            np.savetxt("output/"+out_pre+"_galaxy_dCls_"+inv_nameMap[paramName]+".csv",dCls,delimiter=",")
        else:
            np.savetxt("output/"+out_pre+"_galaxy_dCls_"+paramName+".csv",dCls,delimiter=",")
    print 'End of program'
    '''
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
    '''
if (__name__ == "__main__"):
    main(sys.argv[1:])
