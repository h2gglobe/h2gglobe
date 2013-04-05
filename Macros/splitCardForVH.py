import os
import sys

f = open(sys.argv[1])
y = sys.argv[2]

whFrac2011=[0.65,0.63,0.65,0.63,0.72]
whFrac2012=[0.65,0.63,0.65,0.63,0.64,0.80,0.63,0.83,0.81]
whFrac=[]
if y=="2011":
  whFrac=whFrac2011
elif y=="2012":
  whFrac=whFrac2012
else:
  print 'Unrecognised year', sys.argv[2]

atrel=False
rateLine=False
oldProcesses=['ggH','qqH','VH','ttH','bkg_mass']
newProcesses=['ggH','qqH','WH','ZH','ttH','bkg_mass']
details=[]

for line in f.readlines():
  if line.startswith('shapes VH'):
    print line.replace('VH','WH').replace('wzh','wh'),
    print line.replace('VH','ZH').replace('wzh','zh'),
    continue
  if line.startswith('observation'):
    atrel=True
    print line,
    continue
  if line=='\n':
    print line,
    continue
  if atrel and not line.startswith('process') and not line.startswith('rate'):
    els = line.split()
    print els[0], '   ',
    els = els[1:]
    if els[0]=='lnN': 
      print els[0], '   ',
      els = els[1:]
    for i,el in enumerate(els):
      if (i-2)%len(oldProcesses)==0:
        print el, el,
      else:
        print el,
    print 
    continue
  if atrel and line.startswith('process'):
    if 'VH' in line:
      print line.replace('VH','WH ZH'),
    else:
      print line.replace('-3','-3 -4'),
    continue
  if atrel and line.startswith('rate'):
    els = line.split()
    print els[0], '   ',
    els = els[1:]
    for i, el in enumerate(els):
      if (i-2)%len(oldProcesses)==0:
        print float(el)*whFrac[i//len(oldProcesses)], float(el)*(1.0-whFrac[i//len(oldProcesses)]),
      else:
        print float(el),
    print
    continue
  print line,
     
 

