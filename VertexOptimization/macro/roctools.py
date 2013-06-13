from ROOT import TGraph, TH1F, TSpline3 

class ROCBuilder:
    
    def __init__(self,name,title,sig,bkg):
        self.name = name
        self.title = title
        self.sig = sig
        self.bkg = bkg
        self.buildRoc()

    def buildRoc(self):    
        ## 
        groc = TGraph()
        groc.SetTitle("%s" % self.title )
        groc.GetHistogram().SetTitle("%s; signal eff; background rejection" % self.title )
        groc.SetName("groc_%s" % self.name )
        
        self.roc = TH1F("roc_%s" % self.name, "ROC %s; signal eff; background rejection" % self.title,50,0.,1.)

        nbins = self.sig.GetNbinsX()+1
        sigInt,bkgInt = self.sig.Integral(0,nbins),self.bkg.Integral(0,nbins)
        if sigInt == 0. or bkgInt == 0.:
            return 
        if self.sig.GetMean() > self.bkg.GetMean():
            rng = range(nbins,-1,-1)
        else:
            rng = range(0,nbins+1)
        sigC,bkgC = 0.,0.

        for i in rng:
            sigC += self.sig.GetBinContent(i)
            bkgC += self.bkg.GetBinContent(i)
            
            groc.SetPoint(i, sigC/sigInt, 1. - bkgC/bkgInt )
        
        for i in range( 1, self.roc.GetNbinsX() ):
            x = self.roc.GetBinCenter(i)
            self.roc.SetBinContent( i, groc.Eval(x) )
            
    def getRoc(self):
        return self.roc

    
class ROCIntegrator:
    
    def __init__(self,name,histo):
        self.name = name
        self.histo = histo
        self.rocs = {}
        
    def getRoc(self,order):
        ## print "getRoc %d" % order
        if order == 1:
            return self.histo
        elif order in self.rocs:
            return self.rocs[order]
        else :
            roc = self.getRoc(order-1).Clone("%s_%d" % (self.histo.GetName(), order ))
            roc.Multiply( self.histo )
            self.rocs[order] = roc
            return roc

    def getGraph(self,begin,end):
        gr = TGraph()
        gr.SetName( "eff_%s" %self.name )
        for i in range(begin,end):
            roc = self.getRoc(i)
            gr.SetPoint( gr.GetN(), i, roc.Integral() )
            
        return gr
        
