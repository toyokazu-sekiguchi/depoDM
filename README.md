# Description
Code for calculating energy deposition history due to DM annihilation

# Prerequisites (vesion usued in development)
* Python3 (3.7)
* Pythia8 (8.2.44)
* C++ compiler (gcc Apple LLVM version 10.0.1)

## Python modules
* numpy (1.18.1)
* scipy (1.4.1)
* configparser (4.0.2)
* astropy (4.0)

# Installation
Compilation is required in order to build a python interface for Phythia8 (written in C++)
1. Git clone source file via `git clone https://github.com/toyokazu-sekiguchi/depoDM.git`.
2. Go to `./depoDM/`.
3. Download a copy of Phytia8 from `http://home.thep.lu.se/~torbjorn/pythia8/pythia8244.tgz` and then untar it.
4. Go to `./pythia8244/`, configure with the path to python header `./configure --with-python-include=[path to Python.h]`, and then `make`.
6. Downloads `elec_processed_results.fits` and `phot_processed_results.fits` from Tracy R. Slatyer's webpage `https://faun.rc.fas.harvard.edu/epsilon/detaileddeposition/general/fits/`.
5. Go back to the parent directory. 

# Usage

## Basic usage:
`python3 driver.py params.ini`

This computes the deposition fraction $f_c(z)$ from dm annihilation into five processes, namely 1) hydrogen ionization, 2) helium ionization, 3) Ly-alpha, 4) heating and 5) continuum photons. To get the Energy deposition of each process per time per volume, multiply computed $f_c(z)$ with the rest mass of a dark matter pair (e.g. $2m_{DM}$) times the rate of dark matter annihilation events ($n_{DM} n_{\overline DM}<\sigma v>$). 

### Description of input parameter file
`params.ini` specifies a variety of parameters and consists of five sections:
* [OUTPUT]
  - `root`: Chacters specifying output prefix.
* [COSMOLOGY]
  - `ob`, `odm`, `ode`: Density parameters $\Omega_i h^2$ of baryon, dark matter (assuming no decay) and cosmological constant. Note that actual value of $\Omega_{dm} h^2$ and henceforth $h = \sqrt{\sum_i \Omega_i h^2}$ differ from the input value when decay rate is finite.
  - `nnu`: Effective number of neutrinos. The total number of neutrinos are enhanced by this factor (temperature is fixed to the standard value i.e. $T_\nu = (4/11)^{1/3} T_\gamma$.
  - `mnu`: Sum of neutrino mass in units of eV.
  - `neutrino_hierarchy`: Flag for neutrino mass hierarchy. 1 for normal, 0 for degenerate and -1 for inverted ones.
* [NBODY]
  - `clumpiness`: Path to the text file for look-up table of clumpiness factor; leave it blank if dark matter density is assummed to be uniform.
* [DEPOSITION]
  - `epspath`: Path to the directory where the two `.fits` files from Slatyer's webpage are located.
  - `intype`: Type of injection. 1 for DM annihilation and 2 for DM decay. At present, only annihilation is supported.
* [ANNIHILATION]
  - `mass`: Mass of dark matter in GeV.
  - `mult`: (For future extention) multiplicity factor of DM particle. 1 for Majorana and 2 for Dirac.
  - `mode`: Annihilation products. 1 for two $\gamma$, 2 for $e^+e^-$, 3 for $b\bar{b}$, 4 for $W^+W^-$.

### Role of each python file:
* `const.py`: Definition of units and constants
* `background.py`: Calculation of cosmological background evolution. 
* `inifile.py`: Module for low-level functions used to read `params.ini` file.
* `energy.py`: 
* `ann.py`: 
* `driver.py`: Main function.

## Stage 2: Postprocessing

### Basic usage:
`python3 post.py dist.ini`

This analyses MCMC chain(s) produced in Step 1 and obtain parameter constraints as well as triangle plot of posteior distributions.

### Description of parameter files
`dist.ini` specifies 
* [POSTPROCESS]
  - `postroot`: This specifies output prefix.
  - `paramnames`: Array of parameter names varied in chains. Comma separates items.
  - `paramlabels`: Array of parameter labels in LaTeX format. They are adopted in plotting. Comma separates items.
  - `chains`: Array of chain file(s) to be analysed. Comma separates items.
  - `chainlabels`: Array of chain label(s). They are adopted in plotting. Comma separates items.
  
# Notes
* Flatness is assumed.
* Neutrinos are assumed to consist of three mass eigenstates.
* 4He abundance $Y_p(\omega_b, N_\nu)$ is fitted with a look-up table in `BBN.dat`, which is taken from CLASS, which are originally obtained using the PArthENoPE code (http://parthenope.na.infn.it).
* Recombination history is computed based on HyReC (https://pages.jh.edu/~yalihai1/hyrec/hyrec.html). Hyrec in our code is modified from the original one so that massive neutrinos are incorporated and interface to Python is realized by SWIG.

# Version history
* March 4nd, 2020
  - Restart functionarity is now supported.
  - HyRec wrapper is modified. In the previous version, segmentation faults occur when HyRec is located on a path which includes symbolic links. All the files for look-up tables are now referred by their physical pathes.
* March 2nd, 2020
  - Initial release.

# To-do list 
- [ ] Parameter estimation of derived parameters
