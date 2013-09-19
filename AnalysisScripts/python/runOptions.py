from optparse import OptionParser
parser = OptionParser()
parser.add_option("","--dryRun",action="store_true",dest="dryRun",default=False)
parser.add_option("-i","--inputDat",dest="inputDat")
parser.add_option("-n","--nJobs",dest="nJobs",default=-1)
parser.add_option("-j","--jobId",dest="jobId",default=0)
parser.add_option("-t","--typeRun",type="int",dest="typeRun",default=-1)
parser.add_option("-s","--searchPath",type="string",dest="searchPath",default="common:reduction:baseline:massfac_mva_binned:full_mva_binned")
parser.add_option("-a","--appendSearchPath",dest="postSearchPath",action="append",default=[])
parser.add_option("-p","--prependSearchPath",dest="preSearchPath",action="append",default=[])
parser.add_option("-l","--label",dest="label",default="")
parser.add_option("-v","--verbose",dest="verbose",action='store_true',default=False)
parser.add_option("-w","--watchDutyCycle",dest="watchDutyCycle",action="store_true",default=False)
parser.add_option("--minDutyCycle",dest="minDutyCycle",action="store",type="int",default=0.5)
parser.add_option("--watchDutyCycleAfter",dest="watchDutyCycleAfter",action="store",type="int",default=15)
parser.add_option("--mountEos",dest="mountEos",action="store_true",default=False)


