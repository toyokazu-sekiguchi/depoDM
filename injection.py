import pythia8
import numpy as np
import sys
from scipy import interpolate

# A derived class for (e+ e- ->) GenericResonance -> various final states.
class Sigma1GenRes(pythia8.SigmaProcess):

    def __init__(self): pythia8.SigmaProcess.__init__(self)
    
    # Number of final-state particles.
    def nFinal(self): return 1
    
    # Evaluate sigmaHat(sHat): dummy unit cross section.
    def sigmaHat(self): return 1.
    
    # Select flavour. No colour or anticolour.
    def setIdColAcol(self):
        self.setId(-11,11,999999)
        self.setColAcol(0,0,0,0,0,0)
        
    # Info on the subprocess.
    def name(self): return "GenericResonance"

    def code(self): return 9001

    def inFlux(self): return "ffbarSame"
    
def main07():

    # Pythia generator.
    pythia = pythia8.Pythia()

    # A class to generate the fictitious resonance initial state.
    sigma1GenRes = Sigma1GenRes()

    # Hand pointer to Pythia.
    pythia.setSigmaPtr(sigma1GenRes)

    # Read in the rest of the settings and data from a separate file.
    pythia.readFile("main07.cmnd")

    # Initialization.
    pythia.init()

    # Extract settings to be used in the main program.
    nEvent = pythia.mode("Main:numberOfEvents")
    nAbort = pythia.mode("Main:timesAllowErrors")

    # Histogram particle spectra.
    eGamma = pythia8.Hist("energy spectrum of photons",        100, 0., 250.)
    eE = pythia8.Hist(    "energy spectrum of e+ and e-",      100, 0., 250.)
    eP = pythia8.Hist(    "energy spectrum of p and pbar",     100, 0., 250.)
    eNu = pythia8.Hist(   "energy spectrum of neutrinos",      100, 0., 250.)
    eRest = pythia8.Hist( "energy spectrum of rest particles", 100, 0., 250.)

    # Begin event loop.
    iAbort = 0
    for iEvent in range(0,nEvent):
    
        # Generate events. Quit if many failures.
        if not pythia.next():
            iAbort += 1
            if iAbort < nAbort: continue
            print(" Event generation aborted prematurely, owing to error!")
            break
      
        # Loop over all particles and analyze the final-state ones.
        for i in range(0,pythia.event.size()):
            if pythia.event[i].isFinal():
                idAbs = pythia.event[i].idAbs()
                eI = pythia.event[i].e()
                if idAbs == 22: 
                    eGamma.fill(eI)
                elif idAbs == 11: 
                    eE.fill(eI)
                elif idAbs == 2212:
                    eP.fill(eI)
                elif (idAbs ==12 or idAbs ==14 or idAbs ==16):
                    eNu.fill(eI)
                else:
                    eRest.fill(eI)
                    print(" Error: stable id = %d" % pythia.event[i].id())

    # Final statistics and histograms.
    pythia.stat()
    eGamma *= 2.5 / nEvent
    eE     *= 2.5 / nEvent
    eP     *= 2.5 / nEvent
    eNu    *= 2.5 / nEvent
    eRest  *= 2.5 / nEvent
    #print(eGamma,eE,eP,eNu,eRest)

    # Write Python code that can generate a PDF file with the spectra.
    # Assuming you have Python installed on your platform, do as follows.
    # After the program has run, type "python main07plot.py" (without the " ")
    # in a terminal window, and open "out07plot.pdf" in a PDF viewer.
    hpl = pythia8.HistPlot("main07plot")
    hpl.frame( "out07plot", "Particle energy spectra", "$E$ (GeV)",
               "$(1/N_{\\mathrm{event}}) \\mathrm{d}N / \\mathrm{d}E$ (GeV$^{-1}$)")
    hpl.add( eGamma, "-", "$\\gamma$")
    hpl.add( eE, "-", "$e^{\\pm}$")
    hpl.add( eP, "-", "$p/\\overline{p}$")
    hpl.add( eNu, "-", "$\\nu$")
    hpl.add( eRest, "-", "others")
    # Use logarithmic y scale.
    hpl.plot(True)

if __name__ == '__main__':
    main07()

class Injection:

    def __init__(self,verbose=0):
        self.verbose = verbose

    def SetParams(self,intype,mass,mult,mode,use_prec,nerg,ergmin,ergmax):
        self.intype = intype
        self.mass = mass
        self.eCM = 2*self.mass if self.intype==1 else self.mass
        self.mult = mult
        self.mode = mode
        self.use_prec = use_prec
        self.nerg = nerg
        self.ergmin = ergmin
        self.ergmax = ergmax
        self.dlnerg = np.log(ergmax/ergmin)/self.nerg
        self.erg = np.logspace(np.log10(self.ergmin),np.log10(self.ergmax),self.nerg,base=10)
        if(not use_prec):
            self.eGamma = pythia8.Hist("energy spectrum of photons",        self.nerg, self.ergmin, self.ergmax, True)
            self.eE = pythia8.Hist(    "energy spectrum of e+ and e-",      self.nerg, self.ergmin, self.ergmax, True)
            self.eP = pythia8.Hist(    "energy spectrum of p and pbar",     self.nerg, self.ergmin, self.ergmax, True)
            self.eNu = pythia8.Hist(   "energy spectrum of neutrinos",      self.nerg, self.ergmin, self.ergmax, True)
            self.eRest = pythia8.Hist( "energy spectrum of rest particles", self.nerg, self.ergmin, self.ergmax, True)

    def SetRate(self,rate):
        import const
        if(self.intype==1): # annihilation cross section in m^3/s
            self.sigmav = rate*const.cm**3
        else: # decay rate in 1/s
            self.gamma = rate
            
    def RunPythia(self):
        pythia = pythia8.Pythia()
        sigma1GenRes = Sigma1GenRes()
        pythia.setSigmaPtr(sigma1GenRes)
        pythia.readFile("inj.cmnd")
        s0 = str(self.eCM)
        s1 = str(self.eCM*0.01)
        pythia.readString("Beams:eCM = "+s0) 
        pythia.readString("999999:all = GeneralResonance void 1 0 0 "+s0+" "+s1+" 0. 0. 0.")
        if(self.mode==1):
            pythia.readString("999999:addChannel = 1 1. 101 22 22  !  -> gamma gamma")
        elif(self.mode==2):
            pythia.readString("999999:addChannel = 1 1. 101 11 -11 !  -> e+ e-")
        elif(self.mode==3):
            pythia.readString("999999:addChannel = 1 1. 101 15 -15 !  -> tau+ tau-")
        elif(self.mode==4):
            pythia.readString("999999:addChannel = 1 1. 101 5  -5  !  -> b bbar")
        elif(self.mode==5):
            pythia.readString("999999:addChannel = 1 1. 101 24 -24 !  -> W+ W-")
        elif(self.mode==6):
            pythia.readString("999999:addChannel = 1 1. 101 13 -13 !  -> mu+ mu-")
        else:
            print("error: only gamma, e, tau, b, W and mu are supported")
            sys.exit(1)
        pythia.init()
        
        nEvent = pythia.mode("Main:numberOfEvents")
        nAbort = pythia.mode("Main:timesAllowErrors")
        
        iAbort = 0
        for iEvent in range(0,nEvent):
    
            if not pythia.next():
                iAbort += 1
                if iAbort < nAbort: continue
                print(" Event generation aborted prematurely, owing to error!")
                break
      
            # Loop over all particles and analyze the final-state ones.
            for i in range(0,pythia.event.size()):
                if pythia.event[i].isFinal():
                    idAbs = pythia.event[i].idAbs()
                    eI = pythia.event[i].e()
                    if idAbs == 22: 
                        self.eGamma.fill(eI)
                    elif idAbs == 11: 
                        self.eE.fill(eI)
                    elif idAbs == 2212:
                        self.eP.fill(eI)
                    elif (idAbs ==12 or idAbs ==14 or idAbs ==16):
                        self.eNu.fill(eI)
                    else:
                        self.eRest.fill(eI)
                        print(" Error: stable id = %d" % pythia.event[i].id())
            
        # Final statistics and histograms.
        pythia.stat()
        self.eGamma *=  1/(self.dlnerg*nEvent)
        self.eE     *=  1/(self.dlnerg*nEvent)
        self.eP     *=  1/(self.dlnerg*nEvent)
        self.eNu    *=  1/(self.dlnerg*nEvent)
        self.eRest  *=  1/(self.dlnerg*nEvent)
        if(self.verbose>0): print(self.eGamma,self.eE,self.eP,self.eNu,self.eRest)
            
        if(self.verbose>1):
            # Write Python code that can generate a PDF file with the spectra.
            # Assuming you have Python installed on your platform, do as follows.
            # After the program has run, type "python main07plot.py" (without the " ")
            # in a terminal window, and open "out07plot.pdf" in a PDF viewer.
            hpl = pythia8.HistPlot("main07plot")
            hpl.frame( "out07plot", "Particle energy spectra", "$E$ (GeV)",
                       "$(1/N_{\\mathrm{event}}) \\mathrm{d}N / \\mathrm{d}\ln E$)")
            hpl.add( self.eGamma, "-", "$\\gamma$")
            hpl.add( self.eE, "-", "$e^{\\pm}$")
            hpl.add( self.eP, "-", "$p/\\overline{p}$")
            hpl.add( self.eNu, "-", "$\\nu$")
            #hpl.add( self.eRest, "-", "others")
            # Use logarithmic y scale.
            hpl.plot(False)

    def GetBinnedNumber(self):
        if(self.eCM<10): #this is where Pythia doesn't work
            self.spec_elec = np.zeros(self.nerg) 
            self.spec_phot = np.zeros(self.nerg)
            if(self.mode==1):
                j = int(np.log(self.mass/self.erg[0])/np.log(self.erg[self.nerg-1]/self.erg[0])*(self.nerg-1))
                self.spec_phot[j] = 1.
            elif(self.mode==2):
                j = int(np.log(self.mass/self.erg[0])/np.log(self.erg[self.nerg-1]/self.erg[0])*(self.nerg-1))
                self.spec_elec[j] = 1.
            elif(self.mode==6):
                j = int(np.log(self.mass/3/self.erg[0])/np.log(self.erg[self.nerg-1]/self.erg[0])*(self.nerg-1))
                self.spec_elec[j] = 1.
            else:
                print("error: only mode==1 or 2 is supported at mass<5GeV because of the limitation of Pythia")
                sys.exit(1)

        elif(self.use_prec):
            from addon import function #please set annmode in the form of  'WW', 'bb', 'tautau', 'ee', or '2gam'
            if(self.mode==1):
                annmode = '2gam'
            elif(self.mode==2):
                annmode = 'ee'
            elif(self.mode==3):
                annmode = 'tautau'
            elif(self.mode==4):
                annmode = 'bb'
            elif(self.mode==5):
                annmode = 'WW'
            else:
                print("error: only gamma, e, tau, b and W are supported")
                sys.exit(1)
            s=function.spectra(self.mass,annmode)
            #debug
            #np.savetxt('debug/'+annmode+'0.txt',np.concatenate([s[0][:,None],s[1][:,None],s[2][:,None],s[3][:,None]],axis=1))
            #debug
            #spl_elec = interpolate.make_interp_spline(np.log(s[0]),(s[2]+s[3])*s[0]*self.dlnerg/self.eCM)
            #spl_phot = interpolate.make_interp_spline(np.log(s[0]),s[1]*s[0]*self.dlnerg/self.eCM)
            #self.spec_elec = np.array([spl_elec(np.log(erg)) for erg in self.erg])
            #self.spec_phot = np.array([spl_phot(np.log(erg)) for erg in self.erg])
            spl_elec = interpolate.make_interp_spline(np.log(s[0]),(s[2]+s[3]))
            spl_phot = interpolate.make_interp_spline(np.log(s[0]),s[1])
            self.spec_elec = np.array([spl_elec(np.log(erg)) for erg in self.erg])*self.erg[:]*self.dlnerg/self.eCM
            self.spec_phot = np.array([spl_phot(np.log(erg)) for erg in self.erg])*self.erg[:]*self.dlnerg/self.eCM
            #debug
            #np.savetxt('debug/'+annmode+'1.txt',np.concatenate([self.erg[:,None],self.spec_phot[:,None],self.spec_elec[:,None]],axis=1))
            #debug
            
        else:
            self.RunPythia()
            self.spec_elec = np.array([self.erg[j]*self.eE.getBinContent(j+1) for j in range(self.nerg)])*self.dlnerg/self.eCM
            self.spec_phot = np.array([self.erg[j]*self.eGamma.getBinContent(j+1) for j in range(self.nerg)])*self.dlnerg/self.eCM
