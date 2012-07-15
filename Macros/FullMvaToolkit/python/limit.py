import sys,numpy,os
from optparse import OptionParser

# UserInput
parser=OptionParser()
parser.add_option("-l","--mHMin",default=110, type='float')
parser.add_option("-u","--mHMax",default=150, type='float')
parser.add_option("-s","--mHStep",default=0.5, type='float')
parser.add_option("-r","--runDirectory",default="./mva-datacards-grad/",type='str')
(options,args)=parser.parse_args()

masses = numpy.arange(options.mHMin,options.mHMax+options.mHStep,options.mHStep)

Methods = ["Asymptotic","ProfileLikelihood","MaxLikelihoodFit"]

os.system("mkdir -p %s/Asymptotic %s/ProfileLikelihood %s/MaxLikelihoodFit %s/ExpProfileLikelihood"%(options.runDirectory,options.runDirectory,options.runDirectory,options.runDirectory))

for m in masses:
  print " -------------------------------"
  print " running mass", m
  print " -------------------------------"
  os.system("combine %s/mva-datacard_grad_%3.1f.txt -M Asymptotic --newExpected -m %3.1f"%(options.runDirectory,m,m))
  os.system("combine %s/mva-datacard_grad_%3.1f.txt -M ProfileLikelihood --signif --pvalue -m %3.1f" %(options.runDirectory,m,m))
  os.system("combine %s/mva-datacard_grad_%3.1f.txt -M MaxLikelihoodFit --rMin -10 --rMax 10  -m %3.1f" %(options.runDirectory,m,m))
  for M in Methods: os.system("mv higgsCombineTest.%s.* %s/%s"%(M,options.runDirectory,M))
  os.system("combine %s/mva-datacard_grad_%3.1f.txt -M ProfileLikelihood --signif --pvalue -m %3.1f -t -1 --expectSignal 1" %(options.runDirectory,m,m))
  os.system("mv higgsCombineTest.ProfileLikelihood.* %s/ExpProfileLikelihood"%options.runDirectory)

print "Finished Limits"





