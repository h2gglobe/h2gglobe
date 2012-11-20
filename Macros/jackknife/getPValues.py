#! /usr/bin/python

import sys, os
import ROOT

bestPValues = {}
mus = {}
#workingDir=sys.argv[1]
datacardName = sys.argv[1]
nJobs=sys.argv[2]


for i in xrange(int(nJobs)):
    r = [120+0.5*k for k in range(20)]
    maxMassPValue = (0, 9999.)
    cmd="python PValuesPerCat.py -s part"+str(i)+" -o -p -m 120,120.5,121,121.5,122,122.5,123,123.5,124,124.5,125,125.5,126,126.5,127,127.5,128,128.5,129,129.5,130 -d part"+str(i)+"/"+datacardName+".txt"
    print cmd
    output=os.system(cmd)

    for j in r:
        file = open("part"+str(i)+"/"+datacardName+"/higgsCombinePValue.ProfileLikelihood.mH"+str(j)+".log")
        lines = file.readlines()
        file.close()
    
        for l in lines:
            if ("p-value" in l):
                pvalue = (float(l.split()[-1]))
                if (pvalue < maxMassPValue[1]):
                    maxMassPValue = (j, pvalue)
    bestPValues[i]=maxMassPValue

    # Compute mu
    bestcmd="python PValuesPerCat.py -s part"+str(i)+" -b -m "+str(maxMassPValue[0])+" -d part"+str(i)+"/"+datacardName+".txt"
    print bestcmd
    output=os.system(bestcmd)
    
    thismu={}
    file = open("part"+str(i)+"/"+datacardName+"_BestFit/higgsCombineTest.ChannelCompatibilityCheck.mH"+str(maxMassPValue[0])+".log")
    lines = file.readlines()
    file.close()
    
    for line in lines:
        if "Nominal fit" in line:
            linelist=line.split()
            thismu["best"]=[linelist[5],linelist[6].split("/")]
        elif "Alternate fit" in line:
            linelist=line.split()
            thismu[linelist[-1]]=[linelist[4],linelist[5].split("/")]
    #a = os.popen("combine "+workingDir+"/part"+str(i)+"/"+datacardName+".txt -m "+str(bestPValues[i][0])+" -M ChannelCompatibilityCheck --rMin=-20 --saveFitResult -s -1 -n SignalStrength")
    #lines = a.readlines()
    #for l in lines:
    #    if ("Nominal" in l):
    #        mus[maxMassPValue]=(float(l.split()[5]))
    mus[i]=thismu 
#print mus

mukeys=mus.keys()
mukeys.sort()


file = open(datacardName+"_bestfits.txt","w")
for mukey in mukeys:
    muinfo=mus[mukey]
    
    file.write(str(mukey)+"\t"+str(muinfo["best"][0])+"\t"+str(bestPValues[mukey][0])+"\t"+str(bestPValues[mukey][1])+"\n")
    bestPValues[i]=maxMassPValue

file.close()
            
