import sys
import numpy as np
import ConfigParser
from makeDerivs import getHubbleCosmology,getPowerCamb

def main(argv):
    # Read Config 
    iniFile = "input/testmnu_makeDerivs.ini"
    Config = ConfigParser.SafeConfigParser()
    Config.optionxform = str
    Config.read(iniFile)
    spec = Config.get('general','spec')
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

    # Save starting fiducials
    fidparams = fparams.copy()
    print fidparams

    #testList = Config.get('testmnu','testList').split(',')
    #spread = Config.getfloat('testmnu','spread')
    #testRange = np.arange(1.0-spread/2,1.0+spread/2,spread/6)
    testRange = {}
    testList = []
    for (key, val) in Config.items('testmnu'):
        pmin,pmax,numstep = [float(x) for x in val.split(',')]
        testList.append(key)
        pstep = (pmax-pmin)/numstep
        testRange[key] = np.arange(pmin,pmax,pstep)
        # round the value to about 5 significant digits
        decimal = int(-(np.ceil(np.log10(pmin))-5))
        testRange[key] = [round(x,decimal) for x in testRange[key]]

    print "testRange: ",testRange
    # avoid recalculation with central fiducial set
    #fidDone = 1
    prefix = ''
    for testParam in testList:
        for testVal in testRange[testParam]:
            prefix = testParam+'_'+str(testVal)
            print "Testing ",testParam," = ",testVal

            # Reset to starting fid, then multiply a fraction
            fparams = fidparams.copy()
            fparams[testParam] = testVal

            # Save fiducials 
            print "Calculating and saving fiducial cosmology..."
            if not('H0' in fparams):
                fparams['H0'] = getHubbleCosmology(theta=fparams['theta100'],params=fparams)
            fidCls = getPowerCamb(fparams,spec,AccuracyBoost=AccuracyBoost)
            np.savetxt("output/"+out_pre+spec+"_"+prefix+"_fCls.csv",fidCls,delimiter=",")

            # Calculate and save derivatives
            for paramName in paramList:
                h = stepSizes[paramName]

                print "Calculating forward difference for ", paramName
                pparams = fparams.copy()
                pparams[paramName] = fparams[paramName] + 0.5*h
                if paramName=='theta100':
                    pparams['H0'] = getHubbleCosmology(theta=pparams['theta100'],params=pparams)
                pCls = getPowerCamb(pparams,spec,AccuracyBoost=AccuracyBoost)

                print "Calculating backward difference for ", paramName
                mparams = fparams.copy()
                mparams[paramName] = fparams[paramName] - 0.5*h
                if paramName=='theta100':
                    mparams['H0'] = getHubbleCosmology(theta=mparams['theta100'],params=mparams)
                mCls = getPowerCamb(mparams,spec,AccuracyBoost=AccuracyBoost)

                dCls = (pCls-mCls)/h
                np.savetxt("output/"+out_pre+spec+"_"+prefix+"_dCls_"+paramName+".csv",dCls,delimiter=",")

if (__name__=='__main__'):
    main(sys.argv[1:])
