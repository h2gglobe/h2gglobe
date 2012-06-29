import ROOT
import sys

ROOT.gROOT.ProcessLine( \
   "struct Entry{	\
    double r;		\
    double mh;		\
   };"
)
from ROOT import Entry

def findMaximum(points):
  
  max = -9999.
  maxMass = 0.
  for i in points:
    if (i[0] > max):
      max = i[0]
      maxMass = i[1]

  return maxMass, max
    
def getMass(file):
  try:
   tree = file.Get("limit")
  except:
   return -1
  br = tree.GetBranch("mh")
  c = Entry()
  br.SetAddress(ROOT.AddressOf(c,'mh'))
  tree.GetEntry(0)
  return c.mh


def getOBSERVED(file,entry=0):
  try:
   tree = file.Get("limit")
  except:
   return -1
  br = tree.GetBranch("limit")
  m = tree.GetBranch("mh")
  c = Entry()
  br.SetAddress(ROOT.AddressOf(c,'r'))
  m.SetAddress(ROOT.AddressOf(c,'mh'))
  tree.GetEntry(entry)	
  return c.r,c.mh

files = sys.argv[2:]
files.sort()

gr = ROOT.TGraphAsymmErrors(len(files))

F = [ROOT.TFile(f) for f in files]
p = []
for i,f in enumerate(F):
	P,M = getOBSERVED(f,0)
        p.append((P*P, M))
        #M = getMass(f)

# sort by mass
p = sorted(p, key=lambda i:i[1])

maxMass, max = findMaximum(p)
print max, maxMass
for i,point in enumerate(p):
  gr.SetPoint(i, point[1], max - point[0])
  #print point[1], max - point[0]


func = ROOT.TF1("myfunc", "pol2");
func.SetRange(maxMass *0.99, maxMass *1.01)
gr.Fit("myfunc", "r")
a = gr.GetFunction("myfunc").GetParameter(2)
b = gr.GetFunction("myfunc").GetParameter(1)
vtx = -b/(2*a)
valueAtVtx = gr.GetFunction("myfunc").Eval(vtx)
dX = abs(gr.GetFunction("myfunc").GetX(valueAtVtx + 1.) - vtx)
print "Best Mass: %3.2f +- %3.2f"%(vtx, dX)

out = ROOT.TFile(sys.argv[1],"RECREATE")
out.cd()
gr.SetName("observed")
gr.GetXaxis().SetTitle("m_{H}")
gr.GetYaxis().SetTitle("#Chi^{2}")
l = ROOT.TLine(gr.GetXaxis().GetXmin(), 1, gr.GetXaxis().GetXmax(), 1)
l.Draw("SAME")
#gr.SetGridx(1)
#gr.SetGridy(1)
gr.Write()


	
