#!/bin/env python

from optparse import OptionParser, make_option
import sys, os

def read1D(files, x, igr, title, xtitle):
    ch = ROOT.TChain("limit")
    for f in files:
        print f
        ch.AddFile(f)
    ch.Draw("2*deltaNLL:%s" % (x), "", "" )
    gr = ROOT.gROOT.FindObject("Graph").Clone("gr_%s_%d" % (x,igr) )
    gr.SetTitle(title)
    gr.GetXaxis().SetTitle(xtitle)
    gr.GetYaxis().SetTitle("-2 #Delta LL")

    gr.Sort()

    last = None
    for i in range(gr.GetN(),0,-1):
        if gr.GetX()[i-1] == last:
            gr.RemovePoint(i-1)
        last = gr.GetX()[i-1]
    
    return gr

def main(options, args):

    from rootglobestyle import setTDRStyle
    setTDRStyle()
    
    titles = {
        "RV" : "R_{V}^{#gamma #gamma}"## "( #sigma_{VH} + #sigma_{qqH} ) * BR_{#gamma #gamma} / SM",
        }
    styles = [ (ROOT.kBlue,ROOT.kFullCircle),
               (ROOT.kRed+1,ROOT.kOpenTriangleDown),
               (ROOT.kGreen+2,ROOT.kOpenCircle)
               ]
    objs = []
    graphs = []
    
    for ifile in range(len(options.files)):
        file = options.files[ifile]
        label = "file%d" % ifile
        if len(options.labels)>0:
            label = options.labels.pop(0)

        color,marker = styles.pop(0)
        
        
        if "*" in file:
            files = glob(file)
        else:
            files = [file]

        if len(options.variables) == 1:
            x = options.variables[0]
            gr = read1D( files, x, ifile, label, titles[x] )
            gr.SetLineColor(color)
            gr.SetMarkerColor(color)
            gr.SetMarkerStyle(marker)

            sp = ROOT.GraphToTF1( "mygraph%d" % ifile, gr )
            func = ROOT.TF1("myfunc%d" % ifile,sp,0.,10.,1,"GraphToTF1")
            func.SetParameter(0,0.)
            func.SetLineColor(color)
            gr.GetListOfFunctions().AddLast(func)
            
            
            graphs.append(gr)

            objs.append(sp)
            objs.append(func)
            objs.append(gr)
            
        elif len(options.variables) == 2:
            pass

    
    if len(options.variables) == 1:
        axmin = 999.
        axmax = -999.
        for gr in graphs:
            func = gr.GetListOfFunctions().At(0)
            xmin = func.GetMinimumX()

            eminus = xmin - func.GetX(1.,func.GetXmin(),xmin)
            eplus  = func.GetX(1.,xmin,func.GetXmax()) - xmin

            eminus2 = xmin - func.GetX(2.,func.GetXmin(),xmin)
            eplus2  = func.GetX(2.,xmin,func.GetXmax()) - xmin

            axmin = min(axmin,xmin - eminus2)
            axmax = max(axmax,xmin + eplus2)

            print "%s : %1.3g +%1.2g -%1.2g" % ( gr.GetName(), xmin, eplus , eminus )

        x = options.variables[0]
        canv = ROOT.TCanvas(x,x)
        canv.SetGridx()
        canv.SetGridy()
        leg  = ROOT.TLegend(0.35,0.5,0.7,0.9)
        leg.SetLineColor(ROOT.kWhite)
        leg.SetFillStyle(0)
        objs.append(canv)
        objs.append(leg)
        
        graphs[0].Draw("APL")
        graphs[0].GetXaxis().SetRangeUser(axmin,axmax)
        graphs[0].GetYaxis().SetRangeUser(0.,2.)
        leg.AddEntry(graphs[0],"","pl")
        for gr in graphs[1:]:
            gr.Draw("PL")
            leg.AddEntry(gr,"","pl")
        leg.Draw("same")
        
        for fmt in "C","png","pdf":
            canv.SaveAs( "%s.%s" % ( canv.GetName(), fmt ) ) 
        
    return objs
    
if __name__ == "__main__":
    parser = OptionParser(option_list=[
        make_option("-f", "--file",
                    action="append", dest="files",
                    default=[],
                    ),
        make_option("-l", "--label",
                    action="append", dest="labels",
                    default=[],
                    ),
        make_option("-v", "--variable",
                    action="append", dest="variables",
                    default=[],
                    ),
        ])
    
    (options, args) = parser.parse_args()

    ## sys.argv.append("-b")
    import ROOT

    ROOT.gROOT.LoadMacro("%s/GraphToTF1.C+" % os.path.dirname(sys.argv[0]) )

    objs = main( options, args )

