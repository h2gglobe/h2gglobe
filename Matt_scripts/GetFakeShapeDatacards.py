import os, sys, numpy

outpath=sys.argv[1]
files=sys.argv[2:]
"""
if not os.path.isdir(outpath):
  os.makedirs(outpath)

for f in files:
  text = open(f)
  outtext = open(outpath+'/'+os.path.basename(f),'w')
  for i, line in enumerate(text.readlines()):
    if i==6:
      outtext.write('shapes * * FAKE \n')
      outtext.write(line)
    else:
      outtext.write(line)

for m in numpy.arange(110,150.5,0.5):
  os.system("combine "+outpath+"/mva-datacard_grad_%3.1f.txt -M MaxLikelihoodFit -D data_grad -m %3.1f --saveNorm"%(m,m))
  os.system("mv mlfit.root "+outpath+"/mlfit_%3.1f.root"%m)
os.system("mv higgsCombine* "+outpath)

for m in numpy.arange(110,150.5,0.5):
  os.system("python mlfitNormsToText.py %s/mlfit_%3.1f.root > %s/%3.1fout.txt"%(outpath,m,outpath,m))

"""
for m in numpy.arange(110,150.5,0.5):
  os.system("python plotbymH.py %s FakeShapeResults.root %3.1f"%(outpath,m))

for binM in [110,115,120,125,130,135,140,150]:
  os.system("python plotbybin.py %s FakeShapeResults.root %d"%(outpath,binM)) 
      
