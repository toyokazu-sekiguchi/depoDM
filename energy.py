import astropy.io.fits as iofits
import numpy as np
import background
import sys

class EpsData:

    def __init__(self,f,verbose=0):
        hdul = iofits.open(f)
        if(verbose>0): 
            print("=== header of ",f)
            print(repr(hdul[1].header))

        self.z1out = hdul[1].data["OUTPUT_REDSHIFT"].flatten()
        self.erg = hdul[1].data["ENERGY"].flatten()
        self.erg = np.exp(np.log(10)*self.erg)
        self.z1in = hdul[1].data["INPUT_REDSHIFT"].flatten()
        self.ch = hdul[1].data["CHANNELS"].flatten()
        self.nz1out = len(self.z1out)
        self.nerg = len(self.erg)
        self.nz1in = len(self.z1in)
        self.nch = len(self.ch)
        self.tc = hdul[1].data["DEPOSITION_FRACTIONS_NEW"].reshape(self.nch,self.nz1in,self.nerg,self.nz1out)
        self.fion0 = hdul[1].data['F_ION'].reshape(self.nerg,self.nz1out)
        
        hdul.close()
        if(verbose>0): 
            print("===")

class Deposit:
    
    def __init__(self,verbose=0):
        self.verbose = verbose
        self.species = ["elec","phot"]
        self.nspecies = len(self.species)
        self.epsdata =  []

    def SetParams(self,epspath,mode):
        self.epspath = epspath
        if(not epspath.endswith("/")):
            self.epspath = self.epspath+"/"
        self.mode = mode
        if(self.verbose>0):
            print("\n# energy deposition")
            if(self.mode==1):
                print(" injection type: annihilation")
                self.pow = 6
            elif(self.mode==2):
                print(" injection type: decay")
                self.pow = 3
                print("error: decay is not supported")
                sys.exit(1)
            else:
                print("error: only annihilation is supported")
                sys.exit(1)
        for s in self.species:
            f = self.epspath+s+"_processed_results.fits"
            self.epsdata.append(EpsData(f,verbose=self.verbose))
        self.CheckEpsData()
            
    def CheckEpsData(self):
        self.z1out = self.epsdata[0].z1out
        self.erg = self.epsdata[0].erg
        for n in range(1,self.nspecies):
            if(not((self.epsdata[n].z1out==self.z1out).all() and (self.epsdata[n].erg==self.erg).all())):
                print("error: binning in epsdata disagree")
                sys.exit(1)
        self.nz1out = len(self.z1out)
        self.nerg = len(self.erg)
        self.dlnerg = np.log(self.erg[1]/self.erg[0])
        self.ergmin = self.erg[0]/np.exp(0.5*self.dlnerg)
        self.ergmax = self.erg[self.nerg-1]*np.exp(0.5*self.dlnerg)
        
    def Calcfc(self,dtauda,clumz):
        for n in range(self.nspecies):
            gin = np.array([z1k**(self.pow-5)*dtauda(1/z1k)*clumz(z1k-1) for z1k in self.epsdata[n].z1in])
            gout = np.array([z1i**(self.pow-5)*dtauda(1/z1i) for z1i in self.epsdata[n].z1out])
            self.epsdata[n].fc = np.sum(self.epsdata[n].tc[:,:,:,:]*gin[None,:,None,None]/gout[None,None,None,:],axis=1)

            # check
            #for i in range(self.epsdata[n].nz1out):
            #    for j in range(self.epsdata[n].nerg):
            #        d = np.abs(np.sum(self.epsdata[n].fc[0:2,j,i])/self.epsdata[n].fion0[j,i]-1)
            #        if(d>0.02): print(i,j,d)
                    

