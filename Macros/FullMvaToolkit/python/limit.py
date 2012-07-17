import sys,numpy,os
from optparse import OptionParser

mymasses = []
def cback(option,opt_str,value,parser):
	value = value.split(",")
	for v in value: mymasses.append(float(v))

# UserInput
parser=OptionParser()
parser.add_option("-l","--mHMin",default=110, type='float')
parser.add_option("-u","--mHMax",default=150, type='float')
parser.add_option("-s","--mHStep",default=0.5, type='float')
parser.add_option("-r","--runDirectory",default="./mva-datacards-grad/",type='str')
parser.add_option("","--noas",default=False, action="store_true")
parser.add_option("","--nopl",default=False, action="store_true")
parser.add_option("","--nomf",default=False, action="store_true")
parser.add_option("","--noexppl",default=False, action="store_true")
parser.add_option("","--masses",dest="jobs",action="callback",callback=cback,type='string')

(options,args)=parser.parse_args()

masses = numpy.arange(options.mHMin,options.mHMax+options.mHStep,options.mHStep)
if len(mymasses)>0: masses = mymasses

Methods = []
if not options.noas: Methods.append("Asymptotic")
if not options.nopl: Methods.append("ProfileLikelihood")
if not options.nomf: Methods.append("MaxLikelihoodFit")

os.system("mkdir -p %s/Asymptotic %s/ProfileLikelihood %s/MaxLikelihoodFit %s/ExpProfileLikelihood"%(options.runDirectory,options.runDirectory,options.runDirectory,options.runDirectory))

for m in masses:
  print " -------------------------------"
  print " running mass", m
  print " -------------------------------"
  if not options.noas: 
	os.system("combine %s/mva-datacard_grad_%3.1f.txt -M Asymptotic --newExpected -m %3.1f --rAbsAcc 0.00001 --rRelAcc 0.00001 --minosAlgo=stepping --noFitAsimov"%(options.runDirectory,m,m))
  if not options.nopl: 
	os.system("combine %s/mva-datacard_grad_%3.1f.txt -M ProfileLikelihood --signif --pvalue -m %3.1f" %(options.runDirectory,m,m))
  if not options.nomf: 
	os.system("combine %s/mva-datacard_grad_%3.1f.txt -M MaxLikelihoodFit --rMin -10 --rMax 10  -m %3.1f" %(options.runDirectory,m,m))
  for M in Methods: os.system("mv -v higgsCombineTest.%s.* %s/%s"%(M,options.runDirectory,M))
  if not options.noexppl: 
	os.system("combine %s/mva-datacard_grad_%3.1f.txt -M ProfileLikelihood --signif --pvalue -m %3.1f -t -1 --expectSignal 1" %(options.runDirectory,m,m))
  	os.system("mv -v higgsCombineTest.ProfileLikelihood.* %s/ExpProfileLikelihood"%options.runDirectory)

print "Finished Limits"





