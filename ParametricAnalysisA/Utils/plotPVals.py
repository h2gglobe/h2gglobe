print "Loading Root..."

#import pdb; pdb.set_trace()
from ROOT import *
import os
import sys
import array

gStyle.SetOptStat(000000)
gStyle.SetCanvasBorderMode(0);
gStyle.SetCanvasColor(kWhite);

print "Setting Initial Parameters."
can = TCanvas("Plots","Plots",850,600)
can.SetLogy(True)
can.SetGrid(True)
leg = TLegend(0.6, 0.72, 0.87, 0.87)
leg.SetFillColor(0)
leg.SetBorderSize(1)
mytext = TLatex()
mytext.SetTextSize(0.04)
mytext.SetNDC()
intlumi=4.76

Masses=array.array("f",[x * 0.1 for x in range(1100,1501,5)])
#print Masses
directories = sys.argv[1:-1]

PValues = {}
PValuesError = {}


for directory in directories:
  logfiles = os.popen("/bin/ls "+directory+"/mass-1??*.log").readlines()
  for file in logfiles: 
    filename = file.strip("\n")

    import re
    mo = re.search("-([0-9\.]+)\.log$", filename)
    if not mo:
      print "filename which did not match:",filename
      assert(mo)

    mass = float(mo.group(1))

    value = float(os.popen("grep p-value "+filename+" | awk '{print $3}'").readlines()[0].strip(")\n"))
    error = float(0.0)

    PValues.setdefault(directory, []).append(dict(value = value, error = error, mass = mass))
    # PValuesError.setdefault(directory, {})[mass] = 
    #PValuesError.append(float(os.popen("grep p-value "+filename+" | awk '{print $5}'").readlines()[0].strip(")\n")))

    PValues[directory].sort(key = lambda entry: entry['mass'])

# end of loop over directories

from pprint import pprint
# pprint(PValues)

colors = [ kBlack, kRed, kBlue, kCyan ]
colorIndex = 0



gcSaver = []

printout = True
if printout:
    multigraph = TMultiGraph()
    multigraph.SetTitle(";M_{Higgs} (GeV/c^{2});P-Value Asymptotic")
    #multigraph.SetMinimum(0.005)
    multigraph.SetMinimum(0.0001)
    multigraph.SetMaximum(10)
    
    graphs={}
    for directory in directories:

        Masses            = [ x['mass'] for x in PValues[directory] ]
        thisPValues       = [ x['value'] for x in PValues[directory] ]
        thisPValuesError  = [ x['error'] for x in PValues[directory] ]

        MassError=array.array("f",[0]*len(Masses))
        
        graph=TGraphErrors(len(PValues[directory]),
                           array.array("f",Masses),
                           array.array("f",thisPValues),
                           MassError,
                           array.array("f",thisPValuesError))


        graphs[directory] = graph

        graph.SetLineColor(colors[colorIndex])

        colorIndex = (colorIndex + 1) % len(colors) 
                           

                           
        graph.SetMarkerStyle(21)
        graph.SetLineWidth(2)
        graph.SetFillColor(kWhite)

        # example directory name
        # combine-pvalue-120412-131630-RELOAD/

        mo = re.search("\d+-\d+-(.*)$", directory)
        if not mo:
            name = directory
        else:
            name = mo.group(1)
        
        leg.AddEntry(graph,name,'L')
        multigraph.Add(graph)

    # end of loop over directories

    multigraph.Draw("AC")
    multigraph.GetXaxis().SetRangeUser(110,150)
    leg.Draw("")


    # add sigma lines
    for sigmas in (1,2,3):

      # this gives the convention that 2 sigma = 5% / 2
      y = RooStats.SignificanceToPValue(sigmas)

      line = TLine(110, y, 150, y)
      line.SetLineWidth(4)

      line.Draw()
      gcSaver.append(line)
      
      label = TLatex(110 + 2, y * 1.1, "%d #sigma" % sigmas)
      label.SetTextAlign(11);
      label.Draw()

      gcSaver.append(label)

    mytext.DrawLatex(0.18,0.8,"#splitline{CMS preliminary}{#sqrt{s} = 7 TeV L = %.2f fb^{-1}}"%float(intlumi))


can.SaveAs(sys.argv[-1])
print "Done!"
