import ROOT
import sys,time

NFILES=-1
NIND=-1

plotinfos={}
samples={}

# Set Defaults

xsize=1000
ysize=800
can=ROOT.TCanvas("can","can",10,10,xsize,ysize)

legx1=0.6
legx2=0.95
legy1=0.6
legy2=0.9
legend = ROOT.TLegend(legx1,legy1,legx2,legy2)
dolegend=True

textx1=0.3
textx2=0.6
texty1=0.8
texty2=0.9
plottext = ROOT.TPaveText(textx1,texty1,textx2,texty2,"brNDC");
dotext=True
text=""

Ncol=1
Nrow=1

stackmax=0
stackmin=0
maxscaleup=1.1
linewidth=2
NReBin=1

Debug=False
DebugNew=True

dolog=False
dogridx=False
dogridy=False
dorebin=False
doreplot=True
doxtitle=True
doytitle=True
dooflow=False
douflow=False
doscale=True
docats=True
singalcat=-1


StaticMin=False
StaticMax=False


DoTitles=True



linex=-1



DoBigLegend=False
DoPopSig=True
DoData=True
DoFill=True
DoFillSig=True
DoCats=True
DoStack=True

DoMulti=False
DoWriteAll=False
DoPrintFlows=False
DoPrintBins=False
DoStats=False
DoSumw2=False
DoInt=False
DoRevInt=False




class PlotInfo:
    def __init__(self):
        self.doplot=0
        self.h2d=0
        self.typplot=0
        self.ncat=1
        self.histoncatindtonames=0
        self.nbinsx=1
        self.nbinsy=1
        self.lowlim=0
        self.highlim=0
        self.lowlim2=0
        self.highlim2=0
        self.xaxislabel=""
        self.yaxislabel=""
        self.plotvarname=""
        self.plotvarcatname=""
        self.catid=0
        self.index=-1
        if Debug:
            "New PlotInfo"
        
    def Print(self):
        print "doplot",self.doplot
        print "h2d",self.h2d
        print "typplot",self.typplot
        print "ncat",self.ncat
        print "histoncatindtonames",self.histoncatindtonames
        print "nbinsx",self.nbinsx
        print "nbinsy",self.nbinsy
        print "lowlim",self.lowlim
        print "highlim",self.highlim
        print "lowlim2",self.lowlim2
        print "highlim2",self.highlim2
        print "xaxislabel",self.xaxislabel
        print "yaxislabel",self.yaxislabel
        print "plotvarname",self.plotvarname
        print "plotvarcatname",self.plotvarcatname
        print "catid",self.catid

    def ColsRows(self):
        if self.ncat==1:
            return 1,1
        elif self.ncat==2:
            return 2,1
        elif self.ncat==3:
            return 3,1
        elif self.ncat==4:
            return 2,2
        elif self.ncat==5 or self.ncat==6:
            return 3,2
        elif self.ncat==7 or self.ncat==8:
            return 4,2
        else:
            return 4,5

        

class SampleInfo:
    def __init__(self):
        self.histfilename=""
        self.itype=-99999
        self.inshortnames=""
        #self.infilenames=""
        self.plotsample=1
        self.addtoleg=1
        self.order=-1
        self.color=-1
        self.scale=1.0
        self.displayname=""
        self.configline=""

    def Print(self):
        print "histfilename",  self.histfilename
        print "itype",         self.itype
        print "inshortnames",  self.inshortnames
        print "plotsample",    self.plotsample
        print "addtoleg",      self.addtoleg
        print "order",         self.order
        print "color",         self.color
        print "scale",         self.scale
        print "displayname",   self.displayname
        print "configline",    self.configline
        

    def ParseConfigLines(self):
        if self.configline is not "":
            items=self.configline.split()
            for item in items:
                if item.find("itype=") is not -1:
                    continue
                elif item.find("order=") is not -1:
                    self.order=int(item.split("=")[1])
                elif item.find("plot=") is not -1:
                    self.plotsample=int(item.split("=")[1])
                elif item.find("leg=") is not -1:
                    self.addtoleg=int(item.split("=")[1])
                elif item.find("color=") is not -1:
                    self.color=int(item.split("=")[1])
                elif item.find("marker=") is not -1:
                    self.maker=int(item.split("=")[1])
                elif item.find("scale=") is not -1:
                    self.scale=float(item.split("=")[1])
                elif item.find("displayname=") is not -1:
                    self.displayname=ROOT.TString(item.split("=")[1])
                    self.displayname.ReplaceAll("@"," ")
                else:
                    print "Item",item,"not parsed."
        else:
            self.plotsample=0
            print "No configuration--no plotting"
   


PLOTPROPS=["dolog","dogridx","dogridy","Debug","DebugNew","doreplot","doxtitle","doytitle","dooflow","douflow","dorebin","StaticMin","StaticMax","dolegend","dotext","Debug","DebugNew"]
def FindFunction(option):
    print option
    if option == "START":
        if len(sys.argv) <2:
            print "Please pass configuration"
            sys.exit(0)
        print "Running startup"
        Startup(sys.argv[1])
    elif option.split()[0] == "setleg":
        Setleg(float(option.split()[1]),
        float(option.split()[2]),
        float(option.split()[3]),
        float(option.split()[4]))
    elif option == "ls":
        PrintPlotNames(plotinfos)
    elif option == "samples":
        PrintSamples(samples)
    elif option == "vars":
        PrintAllVariables()
    elif option.split()[0] == "mvlegx":
        MoveLegX(float(option.split()[1]))
    elif option.split()[0] == "mvlegy":
        MoveLegY(float(option.split()[1]))
    elif option.split()[0] == "mvtextx":
        MoveTextX(float(option.split()[1]))
    elif option.split()[0] == "mvtexty":
        MoveTextY(float(option.split()[1]))
    elif option.split()[0] == "singlecat":
        SetCat(int(option.split()[1]))
    elif option.split()[0] == "text":
        NewText(option)
    elif option.split()[0] == "ploton":
        PlotSample(int(option.split()[1]))
    elif option.split()[0] == "j":
        Plot(option.split()[1])
    elif option == "n":
        NextPlot()
    elif option == "p":
        PrevPlot()
    elif option.split()[0]=="rebin":
        SetReBin(option.split()[1])
    elif option.split()[0]=="maxval":
        SetMax(option.split()[1])
    elif option.split()[0]=="minval":
        SetMin(option.split()[1])
    elif option.split()[0]=="resize":
        Resize(option.split()[1],option.split()[2])
    elif option.split()[0] == "write":
        if len(option.split()) == 2:
            Write(str(option.split()[1]))
        else: 
            WriteCat(str(option.split()[1]),int(option.split()[2]))
    elif option in PLOTPROPS:
        Switch(option)
    elif len(option.split())==2:
        itemlist=option.split()
        NewValue(str(itemlist[0]),itemlist[1])
    elif option == "help":
        ListCMDs()
    else:
        print option, "not found"



def ListCMDs():
    print "START:               Read in plotvariables and inputfiles trees"
    print "ls:                  List plots"
    print "samples:             List samples"
    print "j #:                 Standard plot of #th plot in plot list"
    print "p:                   Previous plot"
    print "n:                   Next plot"
    print "rebin #:             Rebin histograms by #"
    print "maxval #:            Set the maximum of the plot (0 turns off)"
    print "minval #:            Set the minimum of the plot (0 turns off)"
    print "ploton itype:        Turns sample with itype on/off"
    print "setleg x1 x2 y1 y2   Moves position of the legend"
    print "variable num:        Change variable to num"
    print "vars:                Prints all global variables"
    print "singlecat #:         Display only # cat (-1 turns off)"
    print "mvlegx/y #:          Move legend by # in x/y"
    print "mvtextx/y #:         Move text by # in x/y"
    print "Switch on/off:      ",PLOTPROPS



### sample editing functions
def PlotSample(itype):
    if samples.has_key(itype):
        if samples[itype].plotsample == 1:
            samples[itype].plotsample=0
        else:
            samples[itype].plotsample=1
        ReDraw()
    else:
        print "no sample",itype
        PrintSamples(samples)
        

def PrintAllVariables():
    keys=globals().keys()
    for var in keys:
        if var == "samples":
            continue
        if var == "plotinfos":
            continue
        print "Name Val",var, globals()[var]


### Plot settings / Canvas Properties
def Switch(option):
    var = globals()[option]
    if Debug:
        print "var",var
    if var == True:
        globals()[option]=False
    else:
        globals()[option]=True
    ReDraw()

        
def NewValue(name,newval):
    if globals().has_key(name):
        globals()[name]=type(globals()[name])(newval)
        ReDraw()
    else:
        print "There is no",name,"to modify."

def SetReBin(nrebin):
    global dorebin
    dorebin = True
    global NReBin
    NReBin = int(nrebin)
    ReDraw()


def SetMin(minval):
    global StaticMin,stackmin
    stackmin=int(minval)
    if stackmin==0:
        StaticMin=False
    else:
        StaticMin=True
    ReDraw()


def SetMax(maxval):
    global StaticMax,stackmax
    stackmax=int(maxval)
    if stackmax==0:
        StaticMax=False
    else:
        StaticMax=True
    ReDraw()


def Resize(x,y):
    global xsize
    xsize = int(x)
    global ysize
    ysize = int(y)
    ReDraw()


def SetLeg(x1,x2,y1,y2):
    global dolegend,legx1,legx2,legy1,legy2
    legx1=x1
    legx2=x2
    legy1=y1
    legy2=y2
    dolegend=True
    ReDraw()


def NewText(fulloption):
    global text
    fulloption=fulloption.lstrip("text")
    text=fulloption.lstrip()
    ReDraw()


#  Information

def PrintPlotNames(plotinfos):
    plotinds = plotinfos.keys()
    plotinds.sort()
    for ind in plotinds:
        print ind, plotinfos[ind].plotvarname

def PrintSamples(samples):
    samplenames = samples.keys()
    for itype in samplenames:
        print "Plotting:",samples[itype].plotsample," itype:",itype,"\t",samples[itype].inshortnames



#  Setup Functions

def Startup(configfilename):
    print "Reading in inputfiles and plotvariables from rootfile"
    print "  and reading plot features from config.dat"
    file=open(configfilename, "r")
    configlines=file.readlines()
    #print configlines
    for line in configlines:
        if Debug:
            print "line",line
        if line.find("file=") is not -1:
            rootline=(line.split("file=")[1]).split()[0]
            global rootfile
            rootfile=ROOT.TFile(rootline)
            break
    print "rootline",rootline
    if rootfile is not "":
        ReadPlotVariables(rootfile)
        ReadInputFiles(rootfile)
        
    else:
        print "need root file in config.dat"
        sys.exit(0)
    
    AssociateConfigToSample(configlines)
    SetupGlobalVariables(configlines)

    for itype in samples:
        samples[itype].ParseConfigLines()


def AssociateConfigToSample(configlines):
    for line in configlines:
        if line.find("name=") is not -1:
            if Debug:
                print "AssociateConfigToSample",line
            items=line.split()
            for item in items:
                if item.find("itype=") is not -1:
                    itype_config=int(item.split("=")[1])
                    if samples.has_key(itype_config):
                        if samples[itype_config].configline=="":
                            samples[itype_config].configline=line
                        elif DebugNew:
                            print "Found another input line for",itype_config
                            print "line",line

    
def SetupGlobalVariables(configlines):
    vars2check=globals().keys()
    for line in configlines:
        if line.find("name=") is not -1:
            continue # this is a files line
        if line.find("text=") is not -1:
            global text
            text=line.lstrip("text=")
            text=text.rstrip("\n")
            print line
            print text
            continue 
        if Debug:
            print "SetupGlobalVariables",line
        for item in line.split():
            if Debug:
                print item
            itemlist=item.split("=")
            if itemlist[0] in vars2check:
                if Debug:
                    print globals()[itemlist[0]]
                    print itemlist,type(itemlist[1])
                if itemlist[1]=='True':
                    globals()[itemlist[0]] = True
                elif itemlist[1]=='False':
                    globals()[itemlist[0]] = False
                else:
                    globals()[itemlist[0]] = type(globals()[itemlist[0]])(itemlist[1])
                if Debug:
                    print globals()[itemlist[0]]


def ReadPlotVariables(rootfile):
    rootfile.cd()
    plotvariables=ROOT.gDirectory.Get("plotvariables")
    entries=plotvariables.GetEntriesFast()
    pvnames=ROOT.TClonesArray("TObjString",entries)
    xaxesnames=ROOT.TClonesArray("TObjString",entries)
    yaxesnames=ROOT.TClonesArray("TObjString",entries)
    
    plotvariables.SetBranchAddress("plotvarnames",ROOT.AddressOf(pvnames))
    plotvariables.SetBranchAddress("xaxislabels",ROOT.AddressOf(xaxesnames))
    plotvariables.SetBranchAddress("yaxislabels",ROOT.AddressOf(yaxesnames))
     
    ientry = plotvariables.LoadTree(0)
    plotvariables.GetEntry(0)
    
    if Debug:
        print "plotvariables.Nvar",plotvariables.Nvar
    for i in range(0,plotvariables.Nvar):
        newplotinfo=PlotInfo()
   
        newplotinfo.index=i 
        newplotinfo.doplot=plotvariables.doplot[i]
        newplotinfo.h2d=plotvariables.h2d[i]
        newplotinfo.typplot=plotvariables.typplot[i]
        newplotinfo.ncat=plotvariables.histoncat[i]
        newplotinfo.histoncatindtonames=plotvariables.histoncatindtonames[i]
        newplotinfo.nbinsx=plotvariables.nbinsx[i]
        newplotinfo.nbinsy=plotvariables.nbinsy[i]
        newplotinfo.lowlim=plotvariables.lowlim[i]
        newplotinfo.highlim=plotvariables.highlim[i]
        newplotinfo.lowlim2=plotvariables.lowlim2[i]
        newplotinfo.highlim2=plotvariables.highlim2[i]
        newplotinfo.xaxislabel=plotvariables.xaxislabels[i].GetString()
        newplotinfo.xaxislabel.ReplaceAll("@"," ")
        newplotinfo.yaxislabel=plotvariables.yaxislabels[i].GetString()
        newplotinfo.yaxislabel.ReplaceAll("@"," ")
        newplotinfo.plotvarname=plotvariables.plotvarnames[i].GetString()
        #newplotinfo.plotvarcatname=
        #newplotinfo.catid=0
        
        plotinfos[i]=newplotinfo
        if Debug:
            plotinfos[i].Print()
        #newplotinfo.Print()


def ReadInputFiles(rootfile):
    rootfile.cd()
    inputfiles=ROOT.gDirectory.Get("inputfiles")
    entries=inputfiles.GetEntriesFast()
    filenames=ROOT.TClonesArray("TObjString",entries)
    shortnames=ROOT.TClonesArray("TObjString",entries)
    inputfiles.SetBranchAddress("infilenames",ROOT.AddressOf(filenames))
    inputfiles.SetBranchAddress("inshortnames",ROOT.AddressOf(shortnames))
    
    ientry = inputfiles.LoadTree(0)
    inputfiles.GetEntry(0)

    NFILES      =inputfiles.nfiles
    NIND        =inputfiles.nindfiles
    intlumi     =inputfiles.intlumi


    for ifile in xrange(NIND):
        #print ifile
        newsample=SampleInfo()
        
        newsample.itype           =inputfiles.itype[ifile]
        newsample.infoind         =inputfiles.infoind[ifile]
        newsample.histoind        =inputfiles.histoind[ifile]
        newsample.inshortnames    =inputfiles.inshortnames[ifile].GetString()
        newsample.infilenames     =inputfiles.infilenames[ifile].GetString()
       
        samples[inputfiles.itype[ifile]]=newsample
 
            

### Ploting Functions 

def Plot(num,printsuffix="",printcat=-1):
    global cur_plot
    cur_plot=plotinfos[int(num)]
    print num, cur_plot.plotvarname
    Ncol,Nrow=cur_plot.ColsRows()

    can.Clear()
    can.SetWindowSize(xsize,ysize)
    can.SetCanvasSize(xsize-4,ysize-28)
    if docats:
        can.Divide(Ncol,Nrow)
   
        for ican in range(1,Ncol*Nrow+1):
            can.cd(ican)
            can.cd(ican).SetLogy(dolog)
            can.cd(ican).SetGridx(dogridx)
            can.cd(ican).SetGridy(dogridy)
    
    
    stacks={}
    stacktypes=["bkg","sig","data"]
    first=1
    stackmaxima={}
    for icat in xrange(cur_plot.ncat):
        stackmaxima[icat]=[]
        for stacktype in stacktypes:
            stacks[stacktype+str(icat)]=MakeStack(stacktype,icat)
            stacks[stacktype+"lines"+str(icat)]=MakeOverlayLines(stacktype, icat)
            stackmaxima[icat].append(stacks[stacktype+str(icat)].GetMaximum())          

    if dolegend:
        SetLegend()

    if dotext:
        SetText()

    if DoTitles:
        can.cd().SetTitle(str(cur_plot.plotvarname))

    if docats:
        cats=xrange(cur_plot.ncat)
    else:
        cats=[]
        cats.append(singlecat)


    for icat in cats:
        stackmaxima[icat].sort()
        
        if StaticMax:
            global stackmax
        else:
            stackmax=stackmaxima[icat][-1]*maxscaleup
        if DebugNew:
            print "stackmaxima",stackmaxima
            print "stackmax",stackmax

        if docats:
            can.cd(icat+1)
      
        stacks["bkg"+str(icat)].Draw("hist")
        stacks["bkg"+str(icat)].SetMaximum(stackmax)
        if StaticMin:
            stacks["bkg"+str(icat)].SetMinimum(stackmin)

        if doxtitle==True:
            stacks["bkg"+str(icat)].GetXaxis().SetTitle(str(cur_plot.xaxislabel))
        else:
            stacks["bkg"+str(icat)].GetXaxis().SetTitle("")
        
        if doytitle==True:
            stacks["bkg"+str(icat)].GetYaxis().SetTitle(str(cur_plot.yaxislabel))
        else:
            stacks["bkg"+str(icat)].GetYaxis().SetTitle("")
            
        
        stacks["bkg"+str(icat)].GetXaxis().SetTitleSize(0.05)
        stacks["bkg"+str(icat)].GetYaxis().SetTitleSize(0.06)
        stacks["bkg"+str(icat)].GetYaxis().SetTitleOffset(0.8)
        #if DoTitles:
        #    stacks["bkg"+str(icat)].SetTitle(str(cur_plot.plotvarname))
        #else:
        stacks["bkg"+str(icat)].SetTitle("")
        
        lineorder = stacks["datalines"+str(icat)].keys()
        lineorder.sort()
        lineorder.reverse()
        for index in lineorder:
            itype=index[1]
            if dolegend and icat==0:
                legend.AddEntry(stacks["datalines"+str(icat)][index],str(samples[itype].displayname),"ep"); 
            #stacks["datalines"+str(icat)][index].Draw("ep same")
        
        lineorder = stacks["siglines"+str(icat)].keys()
        lineorder.sort()
        lineorder.reverse()
        for index in lineorder:
            itype=index[1]
            if dolegend and icat==0:
                legend.AddEntry(stacks["siglines"+str(icat)][index],str(samples[itype].displayname),"l"); 
        
        lineorder = stacks["bkglines"+str(icat)].keys()
        lineorder.sort()
        lineorder.reverse()
        if DebugNew:
            print "lineorder",lineorder
        first=1
        for index in lineorder:
            itype=index[1]
            if DebugNew:
                print "index",index
                print "itype",itype
                print "samples[itype].color",samples[itype].color
            stacks["bkglines"+str(icat)][index].SetLineColor(ROOT.kBlack)
            stacks["bkglines"+str(icat)][index].SetLineWidth(linewidth)
            stacks["bkglines"+str(icat)][index].SetFillStyle(1001)
            stacks["bkglines"+str(icat)][index].SetFillColor(int(samples[itype].color))
            stacks["bkglines"+str(icat)][index].Draw("histsame")
            if dolegend and icat==0:
                legend.AddEntry(stacks["bkglines"+str(icat)][index],str(samples[itype].displayname),"f"); 
        

        stacks["sig"+str(icat)].Draw("histsame")
        stacks["data"+str(icat)].Draw("pesame")
        
        if docats:            
            can.cd(icat+1).Modified()
            can.cd(icat+1).Update()
        else:
            can.cd().Modified()
            can.cd().Update()
            

        if dolegend:
            legend.Draw()

        if dotext:
            plottext.Draw()
    can.cd()    
    can.Update()
    if printsuffix is not "":
        if printcat is not -1:
            can.cd(icat).Print(str(cur_plot.plotvarname)+"_cat"+str(printcat)+"."+printsuffix)
        else:
            can.Print(str(cur_plot.plotvarname)+"."+printsuffix)


def SampleMap(sampletype):
    itypes={}  # key is order, maps to itype
    for sample in samples:
        if samples[sample].plotsample==0:
            continue
        if sampletype == "sig" and sample < 0:
            itypes[samples[sample].order]=samples[sample].itype
        elif sampletype == "data" and sample == 0:
            itypes[samples[sample].order]=samples[sample].itype
        elif sampletype == "bkg" and sample >0:
            itypes[samples[sample].order]=samples[sample].itype

    return itypes


def MakeStack(sampletype, icat):

    itypes=SampleMap(sampletype)
    order = itypes.keys()
    order.sort()

    theStack = ROOT.THStack(sampletype+str(icat),sampletype+str(icat))
    for index in order:
        itype = itypes[index]
        if Debug:
            print "SampleMap  index", index
            print "SampleMap samples[itype].plotsample",samples[itype].plotsample
        if samples[itype].plotsample == 1:
            histname=str(cur_plot.plotvarname)+"_cat"+str(icat)+"_"+str(samples[itype].inshortnames)
            if Debug:
                print "SampleMap histname",histname
            rootfile.cd()
            hist=ROOT.gROOT.FindObject(histname).Clone()
            hist=FormatHist(hist,itype,sampletype)
            theStack.Add(hist)
        
    return theStack


def MakeOverlayLines(sampletype, icat):

    itypes=SampleMap(sampletype)
    if DebugNew:
        print "MakeOverlayLines", itypes
    order = itypes.keys()
    order.sort()
    if Debug:
        print "order",order
   
    # the key is a dictionary from order to itype
    # need itype info for plot properties 
    histmap={}

    first=1
    for index in order:
        if DebugNew:
            print "index",index
        itype = itypes[index]
        if DebugNew:
            print "samples[itype].plotsample",samples[itype].plotsample
            print "samples[itype].addtoleg",samples[itype].addtoleg
        if samples[itype].plotsample == 1:
            histname=str(cur_plot.plotvarname)+"_cat"+str(icat)+"_"+str(samples[itype].inshortnames)
            if DebugNew:
                print "histname",histname
            rootfile.cd()
            if first == 1:
                first=0
                thishist=ROOT.gROOT.FindObject(histname).Clone("stacklines_cat"+str(icat)+"_index"+str(index))
                thishist=FormatHist(thishist,itype,sampletype)
                fullhist=thishist.Clone("fullhist")
            else:
                thishist=ROOT.gROOT.FindObject(histname).Clone("stacklines_cat"+str(icat)+"_index"+str(index))
                thishist=FormatHist(thishist,itype,sampletype)
                thishist.Add(fullhist)
                fullhist.Add(FormatHist(ROOT.gROOT.FindObject(histname),itype,sampletype))
                    
            if  samples[itype].addtoleg is 1:
                if DebugNew:
                    print "thishist.GetMaximum()",thishist.GetMaximum()
                    print "fullhist.GetMaximum()",fullhist.GetMaximum()
                    print "index,itype",index,itype
                # individual histogram properties
                #thishist.SetFillColor(samples[itype].color)
                key=index,itype
                histmap[key]=thishist

    if DebugNew:
        print "len(histmap)",len(histmap)
    return histmap
            

def FormatHist(hist,itype,sampletype):
   
    if doscale: 
        hist.Scale(samples[itype].scale)
    
    if douflow:
        uflow = hist.GetBinContent(0)
        bin1content = hist.GetBinContent(1)
        hist.SetBinContent(1,bin1content+uflow)
    
    if dooflow:
        binN = hist.GetNbinsX()
        oflow = hist.GetBinContent(binN+1)
        binNcontent = hist.GetBinContent(binN)
        hist.SetBinContent(binN,binNcontent+oflow)
    
    if sampletype=="sig":
        hist.SetLineColor(int(samples[itype].color))
        hist.SetLineWidth(linewidth)

    if sampletype=="data":
        hist.SetFillColor(ROOT.kBlack)
        hist.SetFillStyle(0)
        hist.SetMarkerStyle(20)
        hist.SetMarkerSize(1.1)
    
    if dorebin:
        hist.Rebin(NReBin)

    return hist


def SetLegend():
    global legend

    legend = ROOT.TLegend(legx1,legy1,legx2,legy2)
    legend.SetBorderSize(0)
    legend.SetTextFont(62)
    legend.SetLineColor(0)
    legend.SetLineStyle(1)
    legend.SetLineWidth(linewidth)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetFillStyle(0)


def MoveLegX(dx):
    global legx1,legx2
    legx1=legx1+dx
    legx2=legx2+dx
    ReDraw()


def MoveLegY(dy):
    global legy1,legy2
    legy1=legy1+dy
    legy2=legy2+dy
    ReDraw()


def SetText():
    global plottext
   
    plottext = ROOT.TPaveText(textx1,texty1,textx2,texty2,"brNDC");
    plottext.SetTextFont(62);
    plottext.SetTextSize(0.0425)
    plottext.SetBorderSize(0)
    plottext.SetLineColor(0)
    plottext.SetLineStyle(0)
    plottext.SetLineWidth(0)
    plottext.SetFillColor(0)
    plottext.SetFillStyle(0)

    if Debug:
        print text
        print text.find("\\n")
    if text.find("\\n") != -1:
        if Debug:
            print "text",text
        for line in text.split("\\n"):
            if Debug:
                print "line",line
            plottext.AddText(str(line))
    else:
        plottext.AddText(text)


def MoveTextX(dx):
    global textx1,textx2
    textx1=textx1+dx
    textx2=textx2+dx
    ReDraw()


def MoveTextY(dy):
    global texty1,texty2
    texty1=texty1+dy
    texty2=texty2+dy
    ReDraw()

def SetCat(cat):
    global singlecat, docats
    if int(cur_plot.ncat) < int(cat)+1:
        print "Cat",cat,"is larger than cur_plot.ncat",cur_plot.ncat
        print "Setting cat to ncat-1"
        cat=cur_plot.ncat-1
    if cat == -1:
        docats=True
    else:
        docats=False
    singlecat=cat
    ReDraw()


# Plot navigation
    
def PrevPlot():
    if globals().has_key("cur_plot"):
        if DebugNew:
            print "cur_plot.index",cur_plot.index
            print "len(plotinfos.keys())",len(plotinfos.keys())
        if cur_plot.index==0:
            prevnum=len(plotinfos.keys())-1
        else:
            prevnum=cur_plot.index-1
    else:
        prevnum=len(plotinfos.keys())-1
    Plot(prevnum)

def NextPlot():
    if globals().has_key("cur_plot"):
        if DebugNew:
            print "cur_plot.index",cur_plot.index
            print "len(plotinfos.keys())",len(plotinfos.keys())
        if cur_plot.index==len(plotinfos.keys())-1:
            nextnum=0
        else:
            nextnum=cur_plot.index+1
    else:
        nextnum=0
    Plot(nextnum)

def ReDraw():
    if doreplot==True:
        num=cur_plot.index
        Plot(num)

def Write(suffix):
    num=cur_plot.index
    Plot(num,suffix)

def WriteCat(suffix, icat):
    num=cur_plot.index
    Plot(num,suffix,icat)






option="START"
ENDLIST=[".q","q","quit","end","exit"]

while option not in ENDLIST:
    try:
        FindFunction(option)
    except:
        print "Option",option,"failed to run."
    option=str(raw_input("Enter option:  "))


