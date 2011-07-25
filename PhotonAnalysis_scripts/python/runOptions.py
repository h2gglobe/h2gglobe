from optparse import OptionParser
parser = OptionParser()
parser.add_option("","--dryRun",action="store_true",dest="dryRun",default=False)
parser.add_option("-i","--inputDat",dest="inputDat")
parser.add_option("-n","--nJobs",dest="nJobs",default=-1)
parser.add_option("-j","--jobId",dest="jobId",default=0)
