import numpy as np
from scipy import interpolate
import const
import background
import inifile
import ann
import energy
import sys


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
paramsfid = np.array([Ini.ReadFloat(section,"ob"),Ini.ReadFloat(section,"odm"),Ini.ReadFloat(section,"ol"),
                      Ini.ReadFloat(section,"nnu"),Ini.ReadFloat(section,"mnu")])
BG.SetParams(paramsfid)

# clumpiness
section = "NBODY"
path = Ini.ReadString(section,"clumpiness")
print("\n# cluminess")
try:
    tab = np.loadtxt(path)
    clumz = interpolate.make_interp_spline(tab[:,0],tab[:,1])
    print(" cluminess factor is read from ",path)

except:
    print(" table file for cluminess factor is missing; cluminess is ignored:")
    clumz = lambda x:1

# energy deposition
DE = energy.Deposit(verbose=1)
section = "DEPOSITION"
DE.SetParams(Ini.ReadString(section,"epspath"),Ini.ReadInt(section,"intype"))
DE.Calcfc(BG.dtauda,clumz)

# injection spectrum phythia
INJ = ann.Injection(verbose=1)
section = "ANNIHILATION"
INJ.SetParams(Ini.ReadFloat(section,"mass"),Ini.ReadInt(section,"mult"),Ini.ReadInt(section,"mode"),
              DE.nerg,DE.ergmin*const.eV/const.GeV,DE.ergmax*const.eV/const.GeV)
INJ.RunPythia()

spec_elec = np.array([DE.erg[j]*INJ.eE.getBinContent(j+1) for j in range(DE.nerg)])*DE.dlnerg*const.eV/(2*INJ.mass*const.GeV)
spec_phot = np.array([DE.erg[j]*INJ.eGamma.getBinContent(j+1) for j in range(DE.nerg)])*DE.dlnerg*const.eV/(2*INJ.mass*const.GeV)
fz_elec = np.sum(DE.epsdata[0].fc[:,:,:]*spec_elec[None,:,None],axis=1)
fz_phot = np.sum(DE.epsdata[1].fc[:,:,:]*spec_phot[None,:,None],axis=1)
fname = root+"_fz.txt"
np.savetxt(fname,np.transpose(np.concatenate([DE.z1out[None,:],fz_elec[:,:]+fz_phot[:,:]],axis=0)))
#for i in range(DE.nz1out):
#    print(DE.z1out[i],fz_elec[:,i],fz_phot[:,i])
