import numpy as np
from scipy import interpolate
import const
import background
import inifile
import injection
import deposition
import sys
#from matplotlib import pyplot as plt
import therm

args = sys.argv
if(len(args)<2):
    print("error: the number of input parameters is not correct; input command must be")
    print(">$ python driver.py params.ini")
    sys.exit(1)
    
Ini = inifile.IniFile(args[1])
Ini.Dump()

section = "OUTPUT"
root = Ini.ReadString(section,"root")
    
# fiducial cosmology
BG = background.Background(verbose=1)
section = "COSMOLOGY"
paramsfid = np.array([Ini.ReadFloat(section,"ob"),Ini.ReadFloat(section,"odm"),Ini.ReadFloat(section,"ode"),
                      Ini.ReadFloat(section,"nnu"),Ini.ReadFloat(section,"mnu")])
BG.SetParams(paramsfid)

# clumpiness
section = "NBODY"
path = Ini.ReadString(section,"clumpiness")
print("\n# cluminess (ignored when intype==2)")
# clumz is ln(B) as function of ln(1+z) constructed from a lookup table of "z B(z)".
try:
    tab = np.loadtxt(path)
    #clumz = interpolate.make_interp_spline(np.log(1+tab[:,0]),np.log(tab[:,1]))
    clumz = lambda x:0
    print(" cluminess factor is read from ",path)
    
except:
    print(" table file for cluminess factor is missing; cluminess is ignored:")
    clumz = lambda x:0

# energy deposition
DE = deposition.Deposit(verbose=1)
section = "DEPOSITION"
DE.SetParams(Ini.ReadString(section,"epspath"),Ini.ReadInt("INJECTION","intype"))
DE.Calcfc(BG.dtauda,clumz)

# injection spectrum phythia
INJ = injection.Injection(verbose=1)
section = "INJECTION"
mspec = 5
INJ.SetParams(DE.intype,Ini.ReadFloat(section,"mass"),Ini.ReadInt(section,"mult"),Ini.ReadInt(section,"mode"),
              DE.nerg*mspec,DE.erg[0]*const.eV/const.GeV,DE.erg[DE.nerg-1]*const.eV/const.GeV)
INJ.SetRate(Ini.ReadFloat(section,"sigmav") if INJ.intype==1 else Ini.ReadFloat(section,"gamma"))
INJ.GetBinnedNumber()

# IGM evolution
TH = therm.Therm(verbose=1)
TH.SetParams(root)
TH.IntegEnergy(DE,INJ)
TH.ThermInput(BG,DE,INJ)

'''
fz_elec = np.empty([DE.nch,DE.nz1out])
fz_phot = np.empty([DE.nch,DE.nz1out])
for n in range(DE.nch):
    for i in range(DE.nz1out):
        spl = interpolate.interp1d(np.log(DE.erg),DE.epsdata[0].fc[n,:,i],kind='cubic',bounds_error=False,fill_value="extrapolate")
        fc = spl(np.log(INJ.erg*const.GeV/const.eV))
        fz_elec[n,i] = np.sum(fc[:]*INJ.spec_elec[:])
        spl = interpolate.interp1d(np.log(DE.erg),DE.epsdata[1].fc[n,:,i],kind='cubic',bounds_error=False,fill_value="extrapolate")
        fc = spl(np.log(INJ.erg*const.GeV/const.eV))
        fz_phot[n,i] = np.sum(fc[:]*INJ.spec_phot[:])

fname = root+"_fz.txt"
np.savetxt(fname,np.transpose(np.concatenate([DE.z1out[None,:],fz_elec[:,:]+fz_phot[:,:]],axis=0)))
plt.figure()
plt.xlabel("Deposition redshift $(1+z)$")
plt.ylabel("Deposition:injection ratio")
plt.xlim([30,2000])
plt.ylim([0,1.2])
plt.xscale("log")
plt.title("Injection energy = "+str(INJ.mass)+"GeV")
plt.plot(DE.z1out,np.sum(fz_elec[0:5,:]+fz_phot[0:5,:],axis=0),label="Total deposited")
plt.plot(DE.z1out,fz_elec[0,:]+fz_phot[0,:],label="H ionization")
plt.plot(DE.z1out,fz_elec[1,:]+fz_phot[1,:],label="He ionization")
plt.plot(DE.z1out,fz_elec[2,:]+fz_phot[2,:],label="Lyman alpha/Excitation")
plt.plot(DE.z1out,fz_elec[3,:]+fz_phot[3,:],label="Heating")
plt.plot(DE.z1out,fz_elec[4,:]+fz_phot[4,:],label="Continuum")
plt.legend()
plt.savefig(root+"_fz.pdf")

Xion = np.empty(DE.nz1out)
Xexc = np.empty(DE.nz1out)
Xheat = np.empty(DE.nz1out)
for i in range(DE.nz1out):
    a = 1/DE.z1out[i]
    #sigmav = 1e-26 #cm^3/s
    #Gamma = sigmav*BG.odmh2*const.rhoch2/(a*a*a)/(INJ.mass*const.GeV)
    Gamma = 1e-26
    Hinv = a*a*BG.dtauda(a)
    Xion[i] = (fz_elec[0,i]+fz_phot[0,i])*Gamma*Hinv*BG.odmh2/BG.obh2/(1-BG.yp)*const.m_H*const.c**2/const.VH
    Xexc[i] = (fz_elec[2,i]+fz_phot[2,i])*Gamma*Hinv*BG.odmh2/BG.obh2/(1-BG.yp)*const.m_H*const.c**2/const.VH*0.75
    Xheat[i] = (fz_elec[3,i]+fz_phot[3,i])*Gamma*Hinv*BG.odmh2/BG.obh2/(1-BG.yp)*const.m_H*const.c**2/const.eV

BG.UpdateTherm(Xion,Xexc,Xheat)
nz1 = 100
z1start = 2000
z1end = 7
dlnz1 = np.log(z1end/z1start)/(nz1-1)
arr = np.empty([nz1,3])
for i in range(nz1):
    from HyRec import pyrec
    arr[i,0] = z1start*np.exp(dlnz1*i)
    arr[i,1] = pyrec.hyrec_xe(1/arr[i,0])
    arr[i,2] = pyrec.hyrec_tm(1/arr[i,0])
np.savetxt("Tgas.txt",arr)
'''
