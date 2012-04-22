import sys
import os

fname = sys.argv[1:]
path = os.path.dirname(fname[0])
pathout = path+'-novbf'

for file in fname:
  basename = os.path.basename(file)
  print basename, pathout
  f = open(path+'/'+basename)
  out = open(pathout+'/'+basename,'w')
  for line in f.readlines():
    print len(line.split()), line
    if len(line.split())==9:
      for i in range(len(line.split())-1):
        out.write(line.split()[i]+'   ')
    elif len(line.split())==41 or len(line.split())==42:
      for i in range(len(line.split())-5):
        out.write(line.split()[i]+'   ')
    else:
      out.write(line)
    out.write('\n')
      

