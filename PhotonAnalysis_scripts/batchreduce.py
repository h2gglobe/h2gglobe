from optparse import OptionParser
import os
import sys

parser = OptionParser()
parser.add_option("-i", "--inputdirectory", dest="inputdir", help="Location of Data")
parser.add_option("-s", "--outputdirectory", dest="outputdir", help="Final Location of Data")
parser.add_option("-t", "--jobtype", dest="jobtype", help="Type of run: -1 for background MC, 0 for Data, and 1 for Signal")
parser.add_option("-b", "--batchlogdirectory", dest="batchlogdir", default="BatchSubmission_DATE", help="Location of batch logs")
parser.add_option("-q", "--queue", dest="queue", default="1nh", help="Batch queue to submit jobs to")
parser.add_option("--debug", action="store_true", dest="debug", default=False, help="Debug Mode")
parser.add_option("--dryrun", action="store_true", dest="dryrun", default=False, help="Do not submit to batch queue.")
(options, args) = parser.parse_args()

if options.inputdir is None or options.outputdir is None:
	print "inputdirectory and outputdirectory are required."
	sys.exit(1)

def FixDirName(DirName):
	if DirName[len(DirName)-1:len(DirName)]!="/": DirName+="/"
	return DirName

HomeDir = os.popen("pwd").readlines()[0].strip("\n")+"/"
InputDir = FixDirName(options.inputdir)
OutputDir = FixDirName(options.outputdir)

BatchDir=options.batchlogdir
DATE=os.popen("date +%F_%R_%S | tr -s ' ' '_' | tr -s ':' '_'").readlines()[0].strip("\n")
BatchDir=BatchDir.replace("DATE",DATE)
if not os.path.exists(BatchDir): os.popen("mkdir "+BatchDir)

if (options.debug):
	print "Home Directory: "+HomeDir
	print "Input Directory: "+InputDir
	print "Output Directory: "+OutputDir
	print "Batch Directory: "+BatchDir

if os.popen("nsls -d "+OutputDir).close() is not None: os.popen("rfmkdir "+OutputDir)
os.popen("rfchmod 775 "+OutputDir)
InputFiles = os.popen("nsfind "+InputDir+" | grep .root").readlines()

for infile in InputFiles:
	outfile = OutputDir+infile[infile.rfind("/"):].strip("\n")
	logfile = BatchDir+infile.strip("\n")[infile.rfind("/"):]+".log"
	bsubcommand = "bsub -q "+options.queue+" -o "+logfile+" "+HomeDir+"bsubreduce.sh "+infile.strip("\n")+" "+outfile+" "+options.jobtype+" "+HomeDir
	bsubcommand = bsubcommand.replace("//","/")
	if options.debug: print bsubcommand
	if not options.dryrun: os.popen(bsubcommand)

print "Done!"
