import numpy as np
import sys,os
import ConfigParser
import itertools

from cambCall import cambInterface

def getClsCamb(out_pre,spec,templateIni,params,cambRoot=os.environ['AXIONCAMB_DIR'],AccuracyBoost=False,seed=0):
    option = 0 # for regular CAMB call
    CAMB = cambInterface(out_pre,templateIni,cambRoot=cambRoot,option=option,seed=seed)
    for key in params:
        CAMB.setParam(key,params[key])
    CAMB.call(suppress=False)
    if spec == 'lensed':
        Cls = CAMB.getCls(lensed=True)
    else:
        Cls = CAMB.getCls(lensed=False)
    CAMB.done()
    return Cls

def main(argv):
    verbose=True
    nameMap={'H0':'hubble','YHe':'helium_fraction','nnu':'massless_neutrinos','s_pivot':'pivot_scalar','t_pivot':'pivot_tensor','As':'scalar_amp(1)','ns':'scalar_spectral_index(1)','tau':'re_optical_depth','num_massive_neutrinos':'massive_neutrinos','TCMB':'temp_cmb','lmax':'l_max_scalar','kmax':'k_eta_max_scalar'}
    inv_nameMap = {v: k for k, v in nameMap.items()}

    # Read Config
    try:
        iniFile = argv[0]
    except:
        iniFile = "input/makeDerivs_axionCAMB.ini"
    Config = ConfigParser.SafeConfigParser()
    Config.optionxform = str
    Config.read(iniFile)
    
    out_pre = Config.get('general','output_prefix')
    spec = Config.get('general','spec')
    AccuracyBoost = Config.getboolean('general','AccuracyBoost')
    templateIni = Config.get('general','templateIni')
    cambRoot = Config.get('general','cambRoot')
    if not os.path.isfile(templateIni):
        templateIni = cambRoot+Config.get('general','templateIni')
                   
    paramList = []
    fparams = {}
    stepSizes = {}
    fidscript = ''
    
    for (key, val) in Config.items('CAMB'):
        fidscript += key+' = '+val+'\n'
        if key in nameMap:
            key = nameMap[key]
        if ',' in val:
            param, step = val.split(',')
            paramList.append(key)
            fparams[key] = float(param)
            stepSizes[key] = float(step)
        else:
            try:
                fparams[key] = eval(val)
            except:
                fparams[key] = val
            '''                       
            if key in ['l_max_scalar','massive_neutrinos']:
                fparams[key] = int(val)
            else:
                fparams[key] = float(val)               
            '''
    if 'massless_neutrinos' in fparams:
        fparams['massless_neutrinos'] -= fparams['massive_neutrinos']

    # Save fid + stepsize
    fidscript = '[camb]\n'+fidscript
    filename = 'output/'+out_pre+'_'+spec+"_CAMB_fid.csv"
    with open(filename,'w') as tempFile:
        tempFile.write(fidscript)
    
    print fparams
    print paramList,stepSizes
    
    # Save fiducials
    print "Calculating and saving fiducial cosmology..."
    if not('omnuh2' in fparams) and ('mnu' in fparams):
        fparams['omnuh2'] = round(fparams['mnu']/93.14,6)
    fCls = getClsCamb(out_pre,spec,templateIni,fparams,cambRoot=cambRoot,AccuracyBoost=AccuracyBoost)
    np.savetxt("output/"+out_pre+'_'+spec+"_CAMB_fCls.csv",fCls,delimiter=",")

    sys.exit()
    # Calculate and save derivatives
    for paramName in paramList:
        h = stepSizes[paramName]        

        print "Calculating forward difference for ", paramName
        pparams = fparams.copy()
        pparams[paramName] = fparams[paramName] + 0.5*h
        if paramName =='mnu':
            pparams['omnuh2'] = round(pparams['mnu']/93.14,6)
        pCls = getClsCamb(out_pre,spec,templateIni,pparams,cambRoot=cambRoot,AccuracyBoost=AccuracyBoost)
    
    
        print "Calculating backward difference for ", paramName
        mparams = fparams.copy()
        mparams[paramName] = fparams[paramName] - 0.5*h
        if paramName =='mnu':
            mparams['omnuh2'] = round(mparams['mnu']/93.14,6)
        mCls = getClsCamb(out_pre,spec,templateIni,mparams,cambRoot=cambRoot,AccuracyBoost=AccuracyBoost)

        dCls = (pCls-mCls)/h
        if paramName in inv_nameMap:
            np.savetxt("output/"+out_pre+'_'+spec+"_CAMB_dCls_"+inv_nameMap[paramName]+".csv",dCls,delimiter=",")
        else:
            np.savetxt("output/"+out_pre+'_'+spec+"_CAMB_dCls_"+paramName+".csv",dCls,delimiter=",")
    print 'End of program'
if (__name__ == "__main__"):
    main(sys.argv[1:])
