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

'''
Description:
   - First add pyfisher folder to path as $FISHER_DIR


Usage

Note:
--- Be careful to not DOUBLE COUNT when setting fsky (example: if fsky_highEllPlanck = 0.6 - fksy_S4)
'''

def tryLoad(filepath,delimiter=None):
    try:
        return np.loadtxt(filepath,delimiter=delimiter)
    except:
        try:
            return np.loadtxt('output/'+filepath,delimiter=delimiter)
        except:
            return np.loadtxt(os.environ['FISHER_DIR']+'/output/'+filepath,delimiter=delimiter)

class FisherForecast:

    def __init__(self,iniFile,prefix='',dClsRoot={},source='CMB'):

        '''
        Description:
        --- Initialized with an iniFile. See input/fisherDefaults.ini for example.

        Parameters:
        --- iniFile (str) - initialization file (what's in there??) 
        --- prefix  (str) - adiitional identifier that follows derivRoot, that can be specified as a parameter
        --- dClsRoot (set/dictionary) - overwrite the derivRoot for specified parameters

        Examples:
        >>> iniFile = os.environ['FISHER_DIR']+'/input/fisherDebug.ini'
        >>> prefix = ''
        >>> dClsRoot ={'w':'June14_testStepSize_vhhAcc_unlensed_scalar_w_0.04074'}
        >>> F = FisherForecast(iniFile,prefix,dClsRoot)

        
        Return:
        --- FisherForecast self
        '''

        self.Config = ConfigParser.SafeConfigParser()
        #self.Config.optionxform = str
        self.Config.read(iniFile)
        
        self.source = source
        self.prefix = prefix
        self.dClsRoot = dClsRoot
        self.compact = self.Config.getboolean('general','compact')
    
        print iniFile
        self.paramList = self.Config.get('general','paramList').split(',')
        self.numParams = len(self.paramList)
        print self.paramList
        self.totFisher = np.zeros((self.numParams,self.numParams))
        self.globalEllMaxT = self.Config.getint('general','globalEllMaxT')
        self.globalEllMaxP = self.Config.getint('general','globalEllMaxP')
        self.globalEllMaxK = self.Config.getint('general','globalEllMaxK')
        try:
            self.galaxyRoot = self.Config.get('galaxy','galaxyRoot')
            self.galaxy = True
        except:
            print 'galaxyRoot not found!'
            self.galaxy = False
        # Get fiducial setting
        try:
            self.derivRoot = self.Config.get('general','derivRoot')
            if self.prefix != '':
                self.derivRoot += '_'+self.prefix
            fconfig = ConfigParser.SafeConfigParser()
            fconfig.optionxform = str
            try:
                fconfig.read(self.derivRoot+'_fid.csv')
                for (key, val) in fconfig.items('camb'):
                    continue
            except:
                try:
                    fconfig.read('output/'+self.derivRoot+'_fid.csv')
                    fconfig.items('camb')
                except:
                    fconfig.read(os.environ['FISHER_DIR']+'/output/'+self.derivRoot+'_fid.csv')
                    fconfig.items('camb')
            self.fparams={}
            self.stepSizes={}
            for (key, val) in fconfig.items('camb'):
                if ',' in val:
                    param, step = val.split(',')
                    self.fparams[key] = float(param)
                    self.stepSizes[key] = float(step)
                else:
                    self.fparams[key] = float(val)
            print 'Fiducial file found.'
        except:
            print 'No fiducial file.'
            
        # Get Priors
        self.sigPrior = {}
        try:
            for (key,val) in self.Config.items('prior'):
                if key == 'include': continue
                self.sigPrior[key] = float(val)
        except:
            print 'No prior.'
        '''
        try:
            self.sigPrior['tau'] = self.Config.getfloat('general','tauPriorSigma')
        except:
            pass
        '''


    def Noises(self,ell,beamFWHMArcmins,deltaTs,deltaPs,atmT=[0.,0.],atmP=[0.,0.]):
        '''
        A gaussian beam inverted to give you N_ells
        beamFWHMArcmin is beam FWHM in arcminutes
        deltaT and deltaP (rms noise) are in units of microKelvin-arcminutes
        (beam, deltaT, deltaP can be lists for a multi-frequency experiment)
        Atmospheric noise: atmT takes in lknee and alpha for temperature (atmP is for polarization)
        '''
        nTTInv = 0.
        nEEInv = 0.
        atmFactorT = 0.
        atmFactorE = 0.
        
        # atmT consists of lkneeT,alphaT
        if atmT[0] > 1.e-3:
            atmFactorT = (atmT[0]/ell)**(-atmT[1])
        if atmP[0] > 1.e-3:
            atmFactorE = (atmP[0]/ell)**(-atmP[1])
                    
        for beamFWHMArcmin,deltaT,deltaP in zip(beamFWHMArcmins,deltaTs,deltaPs):
            beam = beamFWHMArcmin *np.pi/180./60.
            dT = deltaT *np.pi/180./60.
            dP = deltaP *np.pi/180./60.
            nTT = (dT**2.)*np.exp(ell*(ell+1.)*(beam**2.)/8./np.log(2.))
            nEE = (dP**2.)*np.exp(ell*(ell+1.)*(beam**2.)/8./np.log(2.))
            nTTInv += (1./nTT)
            nEEInv += (1./nEE)
            
        rnTT = (1./nTTInv) * (atmFactorT+1.)
        rnEE = (1./nEEInv) * (atmFactorE+1.)
        
        # make noise very large outside global EllMaxes
        if ell>self.globalEllMaxT: rnTT = 1.e40
        if ell>self.globalEllMaxP: rnEE = 1.e40
            
            
        return rnTT,rnEE
    
    def CovFromVecsSmall(self,Cls,ell,nTT=0.,nEE=0.,nkk=0.,lensing=False,flag=False,galaxy=False,galCls=[],ngalCls=[]):
        '''
        For calculating Fisher through Eq.4 of 1402.4108 ("compact" option)
        Pros: easily extendable to optical cross-correlations
        Cons: doesn't currently reproduce non-compact option exactly
        '''

        TT = Cls[ell,0] + nTT
        EE = Cls[ell,1] + nEE
        TE = Cls[ell,3]
    

        if lensing:
            
            dd = (Cls[ell,4] + nkk)
            dt = (Cls[ell,5])

            #mat = np.array([[dd]])

            
            mat = np.array([[TT,TE,dt],
                           [TE,EE,0.],
                           [dt,0.,dd]])
        else:
            mat = np.array([[TT,TE],
                            [TE,EE]])
            if galaxy:
                #print 'Adding galaxy info'
                #print len(galCls)
                #print ell,ngalCls
                ncol = len(galCls[0,:])
                dim = int((np.sqrt(1.+8.*ncol)-1.)/2.)
                if len(ngalCls) == 0:
                    ngalCls = np.zeros(dim)
                galmat = np.zeros([dim,dim])
                igal = 0
                for i in range(dim):
                    for j in range(dim):
                        if i>j:
                            continue
                        if i == j:
                            galmat[i,j] = galCls[ell,igal] + ngalCls[i]
                        else:
                            galmat[i,j] = galCls[ell,igal]
                            galmat[j,i] = galCls[ell,igal]
                        igal += 1
                mat = block_diag(mat,galmat)
                #mat = block_diag(mat,np.array([galmat[i,j]]))
        # Hard code for just ss
        #mat = np.array([galmat[i,j]])
        #print mat
        return mat

    def CovFromVecs(self,Cls,ell,nTT,nEE,nkk=0.,lensing=False):
        '''
        For calculating Fisher through Eq.A.4 of 1509.0747 (not "compact" option)
        '''
    
        TT = Cls[ell,0] + nTT
        EE = Cls[ell,1] + nEE
        TE = Cls[ell,3]

        f=0.5
        #f=1.0
        if lensing:
            if ell>self.globalEllMaxK or ell<40.:
                nkk = 1.e40
            kk = Cls[ell,4] + nkk
            kt = Cls[ell,5]
            
            mat = np.array([[TT**2.,TE**2.,TT*TE,kt*kt],
                            [TE**2.,EE**2.,EE*TE,0.],
                            [TT*TE,EE*TE,f*(TE**2.+TT*EE),0.],
                            [kt*kt,0.,0.,kk*kk]])
        else:
            mat = np.array([[TT**2.,TE**2.,TT*TE],
                            [TE**2.,EE**2.,EE*TE],
                            [TT*TE,EE*TE,f*(TE**2.+TT*EE)]])
    
        return mat


    def calcFisher(self,verbose=True):

        Config = self.Config
        
        # Loop through each section in config, typically disjoint ell ranges or non-overlapping
        # skies contributing independent Fisher information
        for section in Config.sections():
            #if verbose: print "--- Current Fisher ---",self.totFisher
            #if section=='general': continue
            if not(Config.getboolean(section,'include')):
                #print "Skipping ", section
                continue
            if verbose: print "Calculating fisher for section ", section
            derivRoot = Config.get(section,'derivRoot')
            # Add to derivRoot
            if self.prefix != '':
                derivRoot = derivRoot+'_'+self.prefix
            if section=='BAO':

                self.source += '+BAO'
                sigmafk = [float(x)/1000. for x in Config.get('BAO','sigmafkx1000').split(',')]
        	
                # Calculate Fisher Matrix
                Fisher = np.zeros((len(self.paramList),len(self.paramList)))
                dfks = {}
                for paramName in self.paramList:
                    #dfks[paramName] = np.loadtxt("/astro/u/msyriac/repos/pyfisher/output/"+derivRoot+"_dfk_"+paramName+".csv",delimiter=",")
                    dfks[paramName] = tryLoad(derivRoot+'_dfk_'+paramName+'.csv',',')
                    #dfks[paramName] = np.loadtxt(os.environ['FISHER_DIR']+"/output/"+derivRoot+"_dfk_"+paramName+".csv",delimiter=",")
                zrange = [float(x) for x in Config.get('BAO','redshift').split(',')]
        
                paramCombs = itertools.combinations_with_replacement(self.paramList,2)
                for param1,param2 in paramCombs:
                    if verbose: print "Parameter combination : ", param1,param2
                    i = self.paramList.index(param1)
                    j = self.paramList.index(param2)
                    Fz = 0.		
                    for k in range(0,len(zrange)):
                        dfk1 = dfks[param1][k]
                        dfk2 = dfks[param2][k]
                        Fz += dfk1*dfk2/sigmafk[k]**2.
                    Fisher[i,j] = Fz
                    Fisher[j,i] = Fz
                try:    
                    saveFile = Config.get(section,'saveFisher')
                    np.savetxt(saveFile,Fisher)
                    print "Saved to ", saveFile
                except:
                    pass
                self.totFisher += Fisher
                continue
                
            lmin = Config.getint(section,'lmin')
            lmax = Config.getint(section,'lmax')
            fsky = Config.getfloat(section,'fsky')

            try:
                nFileT = Config.get(section,'noiseFileT')
                ellT, nlTT = np.loadtxt(nFileT,unpack=True)
                fnTT = interp1d(ellT,nlTT,bounds_error=False,fill_value=1.e90)
                nFileP = Config.get(section,'noiseFileP')
                ellP, nlEE = np.loadtxt(nFileP,unpack=True)
                fnEE = interp1d(ellP,nlEE,bounds_error=False,fill_value=1.e90)
                funcNoise = True
                print "Using a noise file"
            except Exception, err:
                traceback.print_exc()
                print "Not using a noise file"
                beamFWHMs = [float(x) for x in Config.get(section,'beamFWHMArcmin').split(',')]
                uKArcminTs = [float(x) for x in Config.get(section,'uKArcminT').split(',')]
                uKArcminPs = [float(x) for x in Config.get(section,'uKArcminP').split(',')]
                try:
                    atmT = [float(x) for x in Config.get(section,'atmT').split(',')]
                    atmP = [float(x) for x in Config.get(section,'atmP').split(',')]
                except:
                    atmT = [0.,0.]
                    atmP = [0.,0.]
                    print "No atmospheric noise found. Set to 0. Add atmT,atmP (read more in Noises func)."
                funcNoise = False
                
            

            

            
            lensing = Config.getboolean(section,'includeLensingAuto')
            nlkkLoc = Config.get(section,'NlkkLocation')
            try:
                galaxy = Config.getboolean(section,'includeGalaxy')
            except:
                print 'Add includeGalaxy=True if want galaxy'
                galaxy = False
            if self.galaxy:
                ngalCls = [float(x) for x in Config.get('galaxy','noise').split(',')]
            else:
                ngalCls = [0]
            fidCls = tryLoad(derivRoot+'_fCls.csv',',')
            #fidCls = np.loadtxt(derivRoot+"_fCls.csv",delimiter=",")

            dCls = {}
            for paramName in self.paramList:
                dfile = derivRoot+'_dCls_'+paramName+'.csv'
                if False: #paramName=='w':
                    dCls['w'] = np.zeros((5477, 6))
                    #dfile = 'July25_highAcc_2pt_szar_step_0.3_unlensed_scalar_dCls_w.csv' # !!!!!
                    # 0 1 3 TT EE TE
                    ssize = "0.01"
                    dCls['w'][:,0] = np.load("deriv_%sstep_tt.npy" %ssize)[:5477]
                    dCls['w'][:,1] = np.load("deriv_%sstep_ee.npy" %ssize)[:5477]
                    dCls['w'][:,3] = np.load("deriv_%sstep_te.npy" %ssize)[:5477]
                else:    
                    dCls[paramName] = tryLoad(dfile,',')
                    #print(dCls[paramName].shape)
                #Cls[paramName] = np.loadtxt(derivRoot+"_dCls_"+paramName+".csv",delimiter=",")
        
            if len(self.dClsRoot) != 0:
                for paramName in self.dClsRoot:
                    dCls[paramName] = tryLoad(self.dClsRoot[paramName]+'_dCls_'+paramName+'.csv',',')
                    #dCls[paramName] = np.loadtxt(self.dClsRoot[paramName]+"_dCls_"+paramName+".csv",delimiter=",")
                    print "Change derivRoot of ",paramName," from ",derivRoot," to ",self.dClsRoot[paramName]

            if False:#galaxy and not(lensing):
                print 'Add Galaxy info.'
                galfidCls = tryLoad(self.galaxyRoot+'_fCls.csv',',')
                #galfidCls = np.loadtxt(self.galaxyRoot+"_fCls.csv",delimiter=",")
                galdCls = {}
                for paramName in self.paramList:
                    galdCls[paramName] = tryLoad(self.galaxyRoot+'_dCls_'+paramName+'.csv',',')
                    #galdCls[paramName] = np.loadtxt(self.galaxyRoot+"_dCls_"+paramName+".csv",delimiter=",")
            else:
                # prevent error
                galfidCls = fidCls
                galdCls = dCls

            ellrange = range(lmin,lmax+1)
            Nlkk = lambda x: 0.
            if lensing or galaxy:
                if verbose: print "lensing = True or galaxy = True, using noise for kk"
                ekk, nl = np.loadtxt(nlkkLoc,unpack=True,usecols=[0,1])
                fnl = interp1d(ekk,nl,bounds_error=False,fill_value=1.e40)
                Nlkk = fnl
        
        
            Cls = []
            nCls = []
            # Loop through each unique parameter combination
            Fisher = np.zeros((self.numParams,self.numParams))
            paramCombs = itertools.combinations_with_replacement(self.paramList,2)
            for param1,param2 in paramCombs:
                #if verbose: print "Parameter combination : ", param1,param2
                i = self.paramList.index(param1)
                j = self.paramList.index(param2)
                Fell = 0.
                for ell in ellrange:

                    
                    if funcNoise:
                        nTT = fnTT(ell)
                        nEE = fnEE(ell)
                    else:
                        nTT,nEE = self.Noises(ell,beamFWHMs,uKArcminTs,uKArcminPs,atmT,atmP)
        
                    if self.compact:

                        # if ell>self.globalEllMaxT or ell>self.globalEllMaxP:
                        #     nTT = 1.e40
                        #     nEE = 1.e40
                        
                        if ell>self.globalEllMaxK or ell<40.:
                            nkk = 1.e40
                        else:
                            nkk = Nlkk(ell)
                        ngalCls[0] = float(nkk)
                        #if ell<40.:
                        #    nkk = 0.
                        #print Nlkk
                        #if (nkk!=0.0) and (nkk!=1.e40):
                            #print nkk,ngalCls
                        Cov = self.CovFromVecsSmall(fidCls,ell,nTT,nEE,nkk=nkk,lensing=lensing,flag=True,galaxy=galaxy,galCls=galfidCls,ngalCls=ngalCls)
                        dCov1 = self.CovFromVecsSmall(dCls[param1],ell,lensing=lensing,galaxy=galaxy,galCls=galdCls[param1])
                        dCov2 = self.CovFromVecsSmall(dCls[param2],ell,lensing=lensing,galaxy=galaxy,galCls=galdCls[param2])
                        if len(Cov) == 1:
                            InvCov = 1./Cov
                            Fell += (2.*ell+1.) * fsky * InvCov*dCov1*InvCov*dCov2 /2.
                        else:
                            InvCov = np.linalg.inv(Cov)
                            Fell += (2.*ell+1.) * fsky * np.trace(np.dot(np.dot(InvCov,dCov1),np.dot(InvCov,dCov2))) /2.
                        
                    else:    
                        Cov = self.CovFromVecs(fidCls,ell,nTT,nEE,nkk=Nlkk(ell),lensing=lensing)
                        if lensing:
                            dCov1 = np.array([dCls[param1][ell,0],dCls[param1][ell,1],
                                        dCls[param1][ell,3],dCls[param1][ell,4]])
                            dCov2 = np.array([dCls[param2][ell,0],dCls[param2][ell,1],
                                            dCls[param2][ell,3],dCls[param2][ell,4]])
                        else:
                            dCov1 = np.array([dCls[param1][ell,0],dCls[param1][ell,1],dCls[param1][ell,3]])
                            dCov2 = np.array([dCls[param2][ell,0],dCls[param2][ell,1],dCls[param2][ell,3]])
                        InvCov = np.nan_to_num(np.linalg.inv(Cov))
                        Fell += (2.*ell+1.) * fsky * np.nan_to_num(dCov1.dot(InvCov).dot(dCov2)) / 2.
        
                    
                Fisher[i,j] = Fell
                Fisher[j,i] = Fell

            try:    
                saveFile = Config.get(section,'saveFisher')
                np.savetxt(saveFile,Fisher)
                print "Saved to ", saveFile
            except:
                pass
            self.totFisher += Fisher
            #print self.totFisher
            
        # Add priors
        FisherPriors = np.zeros((self.numParams,self.numParams))
        for param in self.paramList:
            i = self.paramList.index(param)
            try:
                FisherPriors[i,i] = 1./self.sigPrior[param]**2.
                print '---Including prior for',param,str(self.sigPrior[param])
            except KeyError:
                pass
        if (FisherPriors == 0).all():
            print '---No prior included'
        #print "prior matrix ",FisherPriors
        #print "fisher before adding prior",self.totFisher
        self.totFisher += FisherPriors
        #print "fisher after adding prior",self.totFisher
        #np.savetxt('data/Feb19_DESI_BAO_wcdm.csv',self.totFisher)

        # Final Fisher
        if verbose:    
            print "------- Final Fisher -------\n",self.totFisher
            #print np.linalg.dConfidenceet(self.totFisher)
            try:
                for param in self.paramList:
                    print param, " ," ,self.margSigma(param)
            except:
                pass

        return self.totFisher

    def margSigma(self,param):
        i = self.paramList.index(param)
        return np.sqrt(np.linalg.inv(self.totFisher)[i,i])
        
    def confEllipse(self,param1,param2,confLevel=1,savefig=True,savedata=False,verbose=False):
        script = '[confEllipse]'
        alpha = {1:1.52, 2:2.48, 3:3.41}
        i = self.paramList.index(param1)
        j = self.paramList.index(param2)
        invF = np.linalg.inv(self.totFisher)
        sigmax2 = invF[i,i]
        sigmay2 = invF[j,j]
        sigmaxy = invF[i,j]
        # Ellipse parameters
        if ((sigmay2/sigmax2)<1e-10) or ((sigmax2/sigmay2)<1e-10):
            a = alpha[confLevel]*np.sqrt(max(sigmax2,sigmay2) + sigmaxy**2/max(sigmax2,sigmay2))
            b = alpha[confLevel]*np.sqrt(min(sigmax2,sigmay2) - sigmaxy**2/max(sigmax2,sigmay2))
        else:
            a = alpha[confLevel]*np.sqrt((sigmax2+sigmay2)/2. + np.sqrt( (sigmax2-sigmay2)**2/4. + sigmaxy**2 ))
            b = alpha[confLevel]*np.sqrt((sigmax2+sigmay2)/2. - np.sqrt( (sigmax2-sigmay2)**2/4. + sigmaxy**2 ))
        angle = 1./2.*np.arctan(2.*sigmaxy/(sigmax2-sigmay2))*180./np.pi
        if (sigmax2<sigmay2):
            c = a
            a = b
            b = c
        width = 2.*a
        height = 2.*b
        xcenter = self.fparams[param1]
        ycenter = self.fparams[param2]

        script += '\nsigmax2 = '+str(sigmax2)
        script += '\nsigmay2 = '+str(sigmay2)
        script += '\nsigmaxy = '+str(sigmaxy)     
        script += '\nparams = '+param1+','+param2
        script += '\nconfLevel = '+str(confLevel)
        script += '\nangle = '+str(angle)
        script += '\nxcenter = '+str(xcenter)
        script += '\nycenter = '+str(ycenter)
        script += '\nwidth = '+str(width)
        script += '\nheight = '+str(height)
        if verbose:
            print('Making Confidence Ellipse for '+param1+' and '+param2+' at Confidence Level of '+str(confLevel)+' sigma')
            print(script)
        fileName = 'output/'+self.derivRoot+'_'+self.source+'_confEllipse_'+param1+'_'+param2+'_'+str(confLevel)+'sigma'

        if savedata:
            with open(fileName+'.csv','w') as tempFile:
                tempFile.write(script)
            print('Saved file '+fileName+'.csv')

        if savefig:
            e = Ellipse((xcenter,ycenter),width,height, angle=angle,fill=False)
            fig = plt.figure(figsize=(10,6))
            ax = fig.add_subplot(111)
            ax.add_patch(e)
            ax.plot(xcenter,ycenter,'r*',markersize=16)
            ax.set_xlim([xcenter-width,xcenter+width])
            ax.set_ylim([ycenter-height,ycenter+height])
            ax.set_xticks([xcenter-alpha[confLevel]*np.sqrt(sigmax2),xcenter,xcenter+alpha[confLevel]*np.sqrt(sigmax2)], minor=False)
            ax.set_yticks([ycenter-alpha[confLevel]*np.sqrt(sigmay2),ycenter,ycenter+alpha[confLevel]*np.sqrt(sigmay2)], minor=False)
            ax.set_xlabel(param1)
            ax.set_ylabel(param2)
            #ax.set_xticklabels([)
            plt.grid()
            ax.set_title('Confidence Ellipse for '+param1+' and '+param2+' at Confidence Level of '+str(confLevel)+' sigma')
            plt.savefig(fileName+'.png',format='png')
            print('Saved file '+fileName+'.png')
            plt.close()

def main(argv):

    try:
        iniFile = argv[0]
    except:    
        #iniFile = "input/fisher_s4_temp.ini"
        iniFile = os.environ['FISHER_DIR']+'/input/fisherDebug.ini'
        
    prefix = ''
    #dClsRoot ={'w':'June14_testStepSize_vhhAcc_unlensed_scalar_w_0.04074'}
    dClsRoot ={}
    F = FisherForecast(iniFile,prefix=prefix,dClsRoot=dClsRoot)
    print "Calculating Fisher matrix..."
    FisherMat = F.calcFisher(verbose = True)
    # np.savetxt('data/Feb26_FisherMat_Planck_notau_lens_fsky0.6_lcdm.csv',FisherMat)
    #print "1-sigma error on mnu = " , '{:3.0f}'.format(1.e3*F.margSigma("mnu")) , " meV."
    '''
    for param1 in F.paramList:
        for param2 in F.paramList:
            F.confEllipse(param1,param2,confLevel=1,savefig=True,savedata=False,verbose=True)
    '''
    #F.confEllipse('mnu','tau',confLevel=1,savefig=False,savedata=True,verbose=True)
    #F.confEllipse('omch2','w',confLevel=1,savefig=False,savedata=False,verbose=True)
    #F.confEllipse('mnu','nnu',confLevel=1,savefig=False,savedata=True,verbose=True)
    from orphics import stats
    s4 = stats.FisherMatrix(FisherMat,"H0,ombh2,omch2,tau,As,ns,mnu,w0,wa".split(','),delete_params=['mnu','wa'])
    print(s4.sigmas())


if (__name__ == "__main__"):
    main(sys.argv[1:])
