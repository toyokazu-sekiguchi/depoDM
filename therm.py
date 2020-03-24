import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt
import const

class Therm:

    def __init__(self,verbose=0):
        self.verbose = verbose

    def SetParams(self,root):
        self.root = root
        
    def IntegEnergy(self,DE,INJ):
        self.fz = np.empty([DE.nch,DE.nz1out])
        lnerg = np.log(INJ.erg*const.GeV/const.eV)
        for n in range(DE.nch):
            for i in range(DE.nz1out):
                spl_elec = interpolate.interp1d(np.log(DE.erg),DE.epsdata[0].fc[n,:,i],kind='cubic',bounds_error=False,fill_value="extrapolate")
                spl_phot = interpolate.interp1d(np.log(DE.erg),DE.epsdata[1].fc[n,:,i],kind='cubic',bounds_error=False,fill_value="extrapolate")
                self.fz[n,i] = np.sum(spl_phot(lnerg)*INJ.spec_phot[:]+spl_elec(lnerg)*INJ.spec_elec[:])
        if(self.verbose>0):
            np.savetxt(self.root+"_fz.txt",np.transpose(np.concatenate([DE.z1out[None,:],self.fz[:,:]],axis=0)))
            plt.figure()
            plt.xlabel("Deposition redshift $(1+z)$")
            plt.ylabel("Deposition:injection ratio")
            plt.xlim([30,2000])
            plt.ylim([0,1.2])
            plt.xscale("log")
            plt.title("Injection energy = "+str(INJ.eCM/2)+"GeV")
            plt.plot(DE.z1out,np.sum(self.fz[0:5,:],axis=0),label="Total deposited")
            plt.plot(DE.z1out,self.fz[0,:],label="H ionization")
            plt.plot(DE.z1out,self.fz[1,:],label="He ionization")
            plt.plot(DE.z1out,self.fz[2,:],label="Lyman alpha/Excitation")
            plt.plot(DE.z1out,self.fz[3,:],label="Heating")
            plt.plot(DE.z1out,self.fz[4,:],label="Continuum")
            plt.legend()
            plt.savefig(self.root+"_fz.pdf")

    def ThermInput(self,BG,DE,INJ):
        self.Xion = np.empty(DE.nz1out)
        self.Xexc = np.empty(DE.nz1out)
        self.Xheat = np.empty(DE.nz1out)
        for i in range(DE.nz1out):
            a = 1/DE.z1out[i]
            if(INJ.intype==1):
                sigmav = INJ.sigmav
                Gamma = sigmav*BG.odmh2*const.rhoch2/(a*a*a)/(INJ.mass*const.GeV)/INJ.mult # is this correct?
            elif(INJ.intype==2):
                Gamma = INJ.gamma
            Hinv = a*a*BG.dtauda(a)
            x = Gamma*Hinv*BG.odmh2/BG.obh2/(1-BG.yp)*const.m_H*const.c**2
            self.Xion[i] = self.fz[0,i]*x/const.VH
            self.Xexc[i] = self.fz[2,i]*x/(const.VH*0.75)
            self.Xheat[i] = self.fz[3,i]*x/const.eV
        BG.UpdateTherm(self.Xion,self.Xexc,self.Xheat)
        if(self.verbose>0):
            nz1 = 100
            z1start = 2000
            z1end = 1
            dlnz1 = np.log(z1end/z1start)/(nz1-1)
            arr = np.empty([nz1,3])
            for i in range(nz1):
                from HyRec import pyrec
                arr[i,0] = z1start*np.exp(dlnz1*i)
                arr[i,1] = pyrec.hyrec_xe(1/arr[i,0])
                arr[i,2] = pyrec.hyrec_tm(1/arr[i,0])
            np.savetxt(self.root+"_therm.txt",arr)

