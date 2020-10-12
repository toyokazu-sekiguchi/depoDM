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

    def EvolveTspin(self,BG,DE,INJ):
        # initialize interpolation of collisional coupling rates
        filename = "./data/kappaHH.dat" # based on Table 3 in https://arxiv.org/abs/astro-ph/0608032
        data = np.log(np.loadtxt(filename))
        self.spl_kappaHH = interpolate.make_interp_spline(data[:,0],data[:,1])
        filename = "./data/kappaHe.dat" # based on Table 4 in https://arxiv.org/abs/astro-ph/0608032
        data = np.log(np.loadtxt(filename))
        self.spl_kappaHe = interpolate.make_interp_spline(data[:,0],data[:,1])

        nz1=100
        z1start = 1000
        z1end = 1
        dlnz1 = np.log(z1end/z1start)/(nz1-1)
        arr = np.empty([nz1,5])
        
        nH0 = BG.obh2/(const.m_H*const.c*const.c)*(1-BG.yp)*const.rhoch2
        
        for i in range(nz1):
            from HyRec import pyrec
            z1 = z1start*np.exp(dlnz1*i)
            a = 1/z1
            Tr = const.TCMB/const.kB*z1
            xe = pyrec.hyrec_xe(a)
            Tm = pyrec.hyrec_tm(a)
            Tc = Tm # this should be valid when Ly-alpha is optically thick
            
            nH = nH0*z1*z1*z1
            ne = nH*xe
            nH = nH*(1-xe)
            
            lnTm = np.log(Tm)
            x_c = nH*np.exp(self.spl_kappaHH(lnTm))+ne*np.exp(self.spl_kappaHe(lnTm))
            x_c = x_c*const.E21cm/(const.TCMB*z1)/const.A10
            
            x_alpha = 0
            #x_alpha = (16*np.pi*np.pi*np.pi*const.alphaEM/const.m_e/const.c)*0.4162*4/27/const.A10*const.E21cm/const.TCMB/z1
            #x_alpha = x_alpha* a*a*BG.dtauda(a)/16*(const.c*const.hbar*2*np.pi)/(const.VH*const.VH*9/16)*#(dE/dV/dt)_injection*f_exc
                
            #Tspin = (1+x_c+x_alpha)/(1/Tr+x_c/Tm+x_alpha/Tc)
            Tspin = Tm

            tau_21cm = 3*const.c**3*const.hbar*const.A10*nH/(16*const.kB*const.f21cm**2*Tspin)*BG.dtauda(a)*a**2
                
            arr[i,0] = z1
            arr[i,1] = xe
            arr[i,2] = Tspin
            arr[i,3] = tau_21cm
            #arr[i,3] = 8.6e-3*(1-xe)*(Tr/Tspin)*np.sqrt((BG.obh2+BG.odmh2)/0.15*z1/10)*BG.obh2/0.02
            arr[i,4] = (Tspin-Tr)/z1*tau_21cm
            #arr[i,4] = 23e-3*(1-xe)*(1-Tr/Tspin)*np.sqrt((BG.obh2+BG.odmh2)/0.15*z1/10)*BG.obh2/0.02
        if(self.verbose>0):
            np.savetxt(self.root+"_21cm.txt",arr)
        
        # For EDGES
        z_edges=17
        spl_21cm=interpolate.interp1d(arr[:,0],arr[:,4])
        self.dT21cm_edges = spl_21cm([z_edges])[0]
