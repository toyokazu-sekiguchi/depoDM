# Description
Code for calculating energy deposition history into IGM due to DM annihilation and decay.

# Prerequisites (vesion usued in development)
* Python3 (3.7)
* Pythia8 (8.2.44)
* C++ compiler (gcc Apple LLVM version 10.0.1)
* SWIG (3.0.12)

## Python modules
* numpy (1.18.1)
* scipy (1.4.1)
* configparser (4.0.2)
* astropy (4.0)

# Installation
Compilation is required in order to build a python interface for Phythia8 (written in C++ with python interface based on SWIG).
1. Git clone source file via `git clone https://github.com/toyokazu-sekiguchi/depoDM.git`.
2. Go to `./depoDM/`.
3. Downloads `elec_processed_results.fits` and `phot_processed_results.fits` from Tracy R. Slatyer's webpage `https://faun.rc.fas.harvard.edu/epsilon/detaileddeposition/general/fits/` and then put them in a appropriate directory.  
4. Download a copy of Phytia8 from `http://home.thep.lu.se/~torbjorn/pythia8/pythia8244.tgz` and then untar it.
5. Go to `./pythia8244/`, configure with the path to python header `./configure --with-python-include=[path to Python.h]`, and then `make`.
6. Go back to the parent directory. 
7. Add `./pythia82444/lib/` in `PYTHONPATH`.
8. Go to `./HyRec`.
9. Edit `INC_PY` and `LIB_PY` in `Makefile`, each of which gives the path to `Python.h` or `libpython*.*.so`. Then `make pyrec`.
10. Go back to the parent directory.

# Usage

## Basic usage
`python3 driver.py params.ini`

This computes the deposition fraction $f_c(z)$ from dm annihilation into five processes, namely 1) hydrogen ionization, 2) helium ionization, 3) Ly-alpha, 4) heating and 5) continuum photons. To get the energy deposition of each process per time per volume, multiply computed $f_c(z)$ with the rest mass of a dark matter pair (e.g. $2m_{DM}$) times the rate of dark matter annihilation events ($n_{DM} n_{\overline DM}<\sigma v>$). 

## Description of input parameter file
`params.ini` specifies a variety of parameters and consists of five sections:
* [OUTPUT]
  - `root`: Chacters specifying output prefix.
* [COSMOLOGY]
  - `ob`, `odm`, `ode`: Density parameters $\Omega_i h^2$ of baryon, dark matter (assuming no decay) and cosmological constant. Note that actual value of $\Omega_{dm} h^2$ and henceforth $h = \sqrt{\sum_i \Omega_i h^2}$ differ from the input value when decay rate is finite.
  - `nnu`: Effective number of neutrinos. The total number of neutrinos are enhanced by this factor (temperature is fixed to the standard value i.e. $T_\nu = (4/11)^{1/3} T_\gamma$.
  - `mnu`: Sum of neutrino mass in units of eV.
  - `neutrino_hierarchy`: Flag for neutrino mass hierarchy. 1 for normal, 0 for degenerate and -1 for inverted ones.
* [INJECTION]
  - `mass`: Mass of dark matter in GeV.
  - `intype`: Type of injection. 1 for DM annihilation and 2 for DM decay.
  - `mult`: (For future extention) multiplicity factor of DM particle. 1 for Majorana and 2 for Dirac.
  - `mode`: Annihilation products. 1 for two $\gamma$, 2 for $e^+e^-$, 3 $\tau^+\tau^-$, 4 for $b\bar{b}$, 5 for $W^+W^-$.
  - `sigmav`: Annihilation cross section in cm^3/s. This works only when `intype==1` and is ignored otherwise.
  - `gamma`: Decay rate in 1/s. This works only when `intype==2` and is ignored otherwise.
* [NBODY]
  - `clumpiness`: Path to the text file for look-up table of clumpiness factor; leave it blank if dark matter density is assummed to be uniform.
* [DEPOSITION]
  - `epspath`: Path to the directory where the two `.fits` files from Slatyer's webpage are located.

## Role of each python file
* `const.py`: Definition of units and constants
* `background.py`: Calculation of cosmological background evolution. 
* `inifile.py`: Compilation of low-level functions used to read `.ini` file.
* `deposition.py`: Reading and manipulating the transfer functions of energy deposition from Slatyer's results. 
* `injection.py`: Calculation of spectra ($dN/d\ln E_{kin}$) of photons and electrons/positrons per DM annihilation event based on Pythia8.
* `therm.py`: Calculation of IGM thermal history based on HyRec.
* `driver.py`: Main function.

## Description of Output
* `[root]_fz.txt`: This contains the main results, namely the deposition efficiency $f_c(z)$ for five channels. This text file consists of six columns in the following order:
  - "Redshift" $1+z$
  - Deposition fraction into hydrogen ionization  $f_{H~ion}(z)$
  - Helium ionization efficiency $f_{He~ion}(z)$
  - Ly-alpha $f_{Ly-\alpha}(z)$
  - Heating $f_{heat}(z)$
  - Continuum photons $f_{cont}(z)$.
* `[root]_therm.txt`: This contains the IGM thermal history, consisting of three coluns in the following order:
  - "Redshift" $1+z$
  - Ionization fraction $x_e$.
  - Gas temperature $T_m$.
  
# Notes
* Flat Universe is assumed.
* Neutrinos are assumed to consist of three mass eigenstates.
* 4He abundance $Y_p(\omega_b, N_\nu)$ is fitted with a look-up table in `BBN.dat`, which is taken from CLASS, which are originally obtained using the PArthENoPE code (http://parthenope.na.infn.it).

# Refs
* Deposition transfer function: http://arxiv.org/abs/1506.03812
* Phytia code: http://home.thep.lu.se/Pythia/
* HyRec code: https://pages.jh.edu/~yalihai1/hyrec/hyrec.html

# Version history
* March 24th, 2020
  - DM decay is now supported.
  - IGM evolution (ionization fraction & gas temperature) is computed by integrating modified version of HyRec code. Original HyRec is modified so that python 1) interface is enabled. Due to the application limitation of HyRec, the gas temperature $T_m$ should not be larger than the CMB temperature $T_{\rm CMB}$. This limits the rate of injection, e.g. `sigmav<1e-32` and `gamma 1e-25`. The issue of callback functionarity is circumvented simply by passing arrays.
* March 12th, 2020
  - Initial release.

# To-do (?) list
- [ ] Integration of energy spectra with larger MC samples from Hiroshima san.
- [ ] Extension of recombination calculation into $T_m>T_CMB$, where HyRec halts. In addition, at the moment Peebles's C-factor is not accurately implemented in the coefficient of the excitation contribution in the equation of the ionization fraction.
- [ ] Integration of Python version of 21cmFast. Quickie try failed on my local computer (MacBookPro), probably due to inconsistent setup of gcc. (Default gcc, i.e. clang, is mixed with homebrewed one?)
