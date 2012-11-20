#! /usr/bin/python

import sys

bestPValues = []
workingDir=sys.argv[1]
datacardName = sys.argv[2]
nJobs=sys.argv[3]
for i in xrange(int(nJobs)):
    r = [120+0.5*k for k in range(20)]
    maxMassPValue = (0, 9999.)
    
    for j in r:
        print workingDir+"part"+str(i)+"/"+datacardName+"/higgsCOmbinePValie.ProfileLikelihhod.mH"+str(j)+".log"
        file = open(workingDir+"/part"+str(i)+"/"+datacardName+"/higgsCombinePValue.ProfileLikelihood.mH"+str(j)+".log")
        lines = file.readlines()
        file.close()
    
        for l in lines:
            if ("p-value" in l):
                pvalue = (float(l.split()[-1]))
                print pvalue, maxMassPValue[1]
                if (pvalue < maxMassPValue[1]):
                    maxMassPValue = (j, pvalue)
    bestPValues.append(maxMassPValue)


            
