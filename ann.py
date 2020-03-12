import pythia8
import sys

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
    print(eGamma,eE,eP,eNu,eRest)

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

    def SetParams(self,mass,mult,mode,nerg,ergmin,ergmax):
        self.mass = mass
        self.mult = mult
        self.mode = mode
        self.nerg = nerg
        self.ergmin = ergmin
        self.ergmax = ergmax
        
        self.eGamma = pythia8.Hist("energy spectrum of photons",        self.nerg, self.ergmin, self.ergmax, True)
        self.eE = pythia8.Hist(    "energy spectrum of e+ and e-",      self.nerg, self.ergmin, self.ergmax, True)
        self.eP = pythia8.Hist(    "energy spectrum of p and pbar",     self.nerg, self.ergmin, self.ergmax, True)
        self.eNu = pythia8.Hist(   "energy spectrum of neutrinos",      self.nerg, self.ergmin, self.ergmax, True)
        self.eRest = pythia8.Hist( "energy spectrum of rest particles", self.nerg, self.ergmin, self.ergmax, True)
        
    def RunPythia(self):
        pythia = pythia8.Pythia()
        sigma1GenRes = Sigma1GenRes()
        pythia.setSigmaPtr(sigma1GenRes)
        pythia.readFile("ann.cmnd")
        eCM = str(2*self.mass)
        pythia.readString("Beams:eCM = "+eCM) 
        if(self.mode==1):
            pythia.readString("999999:addChannel = 1 1. 101 22 22  !  -> gamma gamma")
        elif(self.mode==2):
            pythia.readString("999999:addChannel = 1 1. 101 11 -11 !  -> e+ e-")
        elif(self.mode==3):
            pythia.readString("999999:addChannel = 1 1. 101 5  -5  !  -> b bbar")
        elif(self.mode==4):
            pythia.readString("999999:addChannel = 1 1. 101 24 -24 !  -> W+ W-")
        else:
            print("error: only gamma, e, b and W are supported")
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
        self.eGamma *= 2.5 / nEvent
        self.eE     *= 2.5 / nEvent
        self.eP     *= 2.5 / nEvent
        self.eNu    *= 2.5 / nEvent
        self.eRest  *= 2.5 / nEvent
        if(self.verbose>0): print(self.eGamma,self.eE,self.eP,self.eNu,self.eRest)
            
            
