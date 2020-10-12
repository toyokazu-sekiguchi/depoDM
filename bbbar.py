import numpy as np
from scipy import optimize
import driver
import os

debug = False

arr_dT21cm = np.array([-0.05,-0.1,-0.15,-0.2])
num_dT21cm = len(arr_dT21cm)

#mass_low = 5.7e0
#mass_high = 6.2e3
mass_low = 20.e0
mass_high = 6.2e3

sigmavm_low = 1e-30

num_mass = 20
dln_sigmav = 0.7

arr_mass = np.logspace(np.log10(mass_low),np.log10(mass_high),num=num_mass,base=10)
arr_sigmav = np.empty(num_dT21cm)

f0 = "constraints/bbbar_default_params.ini"
with open(f0) as f:
    contents0 = f.readlines()

def run_driver(mass,sigmav):
    str_mass = str(mass)
    str_sigmav = str(sigmav)
        
    f1 = "constraints/bbbar/bbbar_"+str_mass+"GeV_"+str_sigmav+"_params.ini"
    contents1=[]
    for s in contents0:
        contents1.append(s.replace("XXX_YYY",str_mass+"GeV_"+str_sigmav+"cm^3sec^{-1}").replace("XXX",str_mass).replace("YYY",str_sigmav))
                
    with open(f1,"w") as f:
        f.writelines(contents1)
        
    return driver.main(f1,0)
    

for mass in arr_mass:

    sigmav_low = sigmavm_low*mass
    sigmav = sigmav_low
    while True:
        
        sigmav *= np.exp(dln_sigmav)
        sigmav_high = sigmav
        dT21cm_edges = run_driver(mass,sigmav)
        if(debug): print(mass,sigmav,dT21cm_edges)
        if(dT21cm_edges>arr_dT21cm[0]):
            break

    arr_sigmav[:] = 0
        
    for i in range(num_dT21cm):
        def func(ln_sigmav,mass,dT21cm):
            return run_driver(mass,np.exp(ln_sigmav))-dT21cm
        
        result = optimize.root_scalar(func,args=(mass,arr_dT21cm[i]),bracket=(np.log(sigmav_low),np.log(sigmav_high)))
        arr_sigmav[i] = np.exp(result.root)

    print(mass,arr_sigmav[0],arr_sigmav[1],arr_sigmav[2],arr_sigmav[3])
    
