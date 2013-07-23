#!/usr/bin/env python
import sys, os, commands,subprocess

class Conf:
	def __init__(self):
		pass

	def __str__(self):
		rt = ""
		for k,v in self.__dict__.iteritems():
			rt += "%s = %s\n" % (k,v)
		return rt
	
	def set_macro(self,var,val):
		self.expdict_[var]=val % self.expdict_

	def read_histfile(self,line):
		self.histline = line
		split_line = line.split()
		for sp in split_line:
			val = sp.split("=")
			
			if val[0] == "output":
				self.outfile = tuple(str(val[1]).rsplit(".",1))
			
			elif val[0] == 'histfile':
				histfile, extension = str(val[1]).rsplit(".",1)
				self.histfile = os.path.basename(histfile), extension
				self.histdir  = os.path.dirname(histfile)
		


if __name__  == "__main__":
	from optparse import OptionParser
	parser = OptionParser()
	parser.add_option("-i","--inputDat",dest="inputDat",help="default: %default")
	parser.add_option("-n","--nJobs",type="int", dest="nJobs",default=-1,help="default: %default")
	parser.add_option("-o","--outputScript",dest="outputScript",default="",help="default: %default")
	parser.add_option("-l","--label",dest="label",default="",help="default: %default")
	parser.add_option("","--runIC",dest="runIC",default=False, action="store_true",help="default: %default")
	parser.add_option("-u","--user",dest="user",default="",help="default: %default")
	parser.add_option("-a","--addfile",dest="addfiles",action="append",default=[],help="default: %default")
	parser.add_option("-N","--notgz",dest="notgz",action="store_true",default=False,help="default: %default")
	parser.add_option("","--combine",dest="combine",action="store_true",default=False, help="make combiner.py jobs default: %default")
	parser.add_option("--skipData",dest="skipData",action="store_true",default=False,help="default: %default")
	parser.add_option("--skipSig",dest="skipSig",action="store_true",default=False,help="default: %default")
	parser.add_option("--skipBkg",dest="skipBkg",action="store_true",default=False,help="default: %default")
	parser.add_option("--onlyData",dest="onlyData",action="store_true",default=False,help="default: %default")
	parser.add_option("--onlySig",dest="onlySig",action="store_true",default=False,help="default: %default")
	parser.add_option("--onlyBkg",dest="onlyBkg",action="store_true",default=False,help="default: %default")
	parser.add_option("-w","--watchDutyCycle",dest="watchDutyCycle",action="store_true",default=False,help="default: %default")
	parser.add_option("-m","--mountEos",dest="mountEos",action="store_true",default=False,help="default: %default")

	(options,args)=parser.parse_args()

	if options.onlyData:
		options.skipData = False
		options.skipSig = True
		options.skipBkg = True
		options.label += "_data"
	elif options.onlySig:
		options.skipData = True
		options.skipSig = False
		options.skipBkg = True
		options.label += "_sig"
	elif options.onlyBkg:
		options.skipData = True
		options.skipSig = True
		options.skipBkg = False
		options.label += "_bkg"

        if options.combine and not options.runIC : sys.exit("Sorry, no support for submitting combiner jobs at CERN available yet!")
	if options.outputScript == "":
		options.outputScript = "%s/sub" % options.label
	if not options.outputScript.endswith("sub"):
		options.outputScript = ("%s/sub" % options.outputScript).replace("//","/")
	
	# Check IC user configs:
	if options.runIC:
		print "Running at IC-HEP"
		if options.user=="" : 
			sys.exit("If you are an IC user, you need to defined options -u username")
	# In case these guys are going to castor, writing out a handy comining file
	p = open(options.inputDat,"r")
	cfg = Conf() 
	datfile=""
	
	files = [ "" for i in xrange(options.nJobs)]
	has_polder = False
	in_comment = False
	for line in p.read().split("\n"):
		if line.startswith("->"):
			in_comment = not in_comment

		if in_comment or line.startswith("#"):
			datfile += "%s\n" % line
			continue

		line = line.replace("%(label)s",options.label)

		if "histfile" in line:
			cfg.read_histfile(line)
			if( cfg.histdir != "" ):
				line = line.replace(cfg.histdir,"$histdir").replace("$histdir/","$histdir")

		keep = False
		if "typ" in line:
			for sp in line.split(" "):
				if "typ=" in sp:
					typ = int(sp.split("=")[-1])
					if typ == 0 and not options.skipData or typ < 0 and not options.skipSig or typ > 0 and not options.skipBkg:
						keep = True
		if line.startswith("split"):
			if keep:
				files.append("%s\n" % line.replace("split ",""))
			if not has_polder:
				datfile += "$input_files"
				has_polder = True
		elif "typ" in line:
			if keep:
				for i in xrange(options.nJobs):
					files[i] += "%s\n" % line
			if not has_polder:
				datfile += "$input_files"
				has_polder = True
		else:
			datfile += "%s\n" % line

	mydir=os.getcwd()
	scriptdir=os.path.dirname(options.outputScript)
	os.system("mkdir -p %s" % scriptdir)

	
	domain = commands.getoutput("hostname -d")
	knownDomains = { "cern.ch" : "root://eoscms//eos/cms" }
	knownDomain = domain in knownDomains
	if domain == "cern.ch":
		atCern = True
	mkdir="mkdir"
	cp="cp -pv"
	prependToStore=""
	if not options.runIC:
	  if cfg.histdir.startswith("/castor"):
		mkdir="rfmkdir"
		cp="rfcp"
	  elif cfg.histdir.startswith("/store"):
		cp="cmsStage -f"
		mkdir="cmsMkdir"
		if knownDomain:
			prependToStore = knownDomains[domain]
	if cfg.histdir=="":
		if os.path.isabs(scriptdir):
			cfg.histdir=scriptdir
		else:
			cfg.histdir=os.path.join(os.getcwd(),scriptdir)

	outnam=os.path.join(scriptdir, "filestocombine_%s" % os.path.basename(options.inputDat))
	g = open(outnam,"w+") 
	if options.runIC:
	  filestocmb = ""
	  for i in xrange(len(files)):
                if cfg.histdir=='':cfg.histdir='./'
		filestocmb += "typ=99999 Fil=%s/%s_%d.%s\n"  %(cfg.histdir,cfg.histfile[0],i,cfg.histfile[1])
	  g.write( datfile.replace("$input_files",filestocmb).replace("$histdir",scriptdir ) )
	  g.close()	
		
	if not options.runIC:
	  filestocmb = ""
	  for i in xrange(len(files)):
	        if knownDomain:
			if prependToStore != "":
				fil = "%s/%s_%d.%s" % ( prependToStore, os.path.join(cfg.histdir,cfg.histfile[0]), i, cfg.histfile[1] )
			else:
				fil = "%s_%d.%s" % ( os.path.join(cfg.histdir,cfg.histfile[0]), i, cfg.histfile[1] )
		else:
			fil = commands.getoutput("cmsPfn %s_%d.%s" % ( os.path.join(cfg.histdir,cfg.histfile[0]), i, cfg.histfile[1] ))
		filestocmb += "typ=99999 Fil=%s\n" % fil
	  g.write( datfile.replace("$input_files",filestocmb).replace("$histdir", "%s/" % scriptdir ) )
	  g.close()	

	tmpnam = os.path.join(os.path.dirname(options.inputDat), "tmp_%s" % os.path.basename(options.inputDat))
	tmp  = open("%s" % tmpnam, "w+")
	idat = open(options.inputDat, "r")
	if options.runIC:tmp.write(idat.read())
	else:	tmp.write( idat.read() )
	

	idat.close()
	tmp.close()
	if not os.path.isfile("%s.pevents" % tmpnam) and not options.combine:
		print "Generating the pevents file...",
		print "python fitter.py -l %s -i %s --dryRun >& %s.log\n" % (options.label,tmpnam,tmpnam)
		os.system("python fitter.py -l %s -i %s --dryRun >& %s.log" % (options.label,tmpnam,tmpnam) )
		print "Done. Check %s.log for errors" % tmpnam

	os.system("%s %s" % ( mkdir, cfg.histdir) )
	if not options.notgz:
		filedir = os.path.dirname(options.inputDat)
		searchpath = "common reduction baseline massfac_mva_binned full_mva_binned jetanalysis photonjet"
		os.system("tar zcf %s.tgz *.py $(find %s -name \*.dat -or -name \*.py) aux common python %s" %
			  (options.outputScript, searchpath, filedir ) )
		os.system("%s %s.tgz %s" % ( cp,  options.outputScript, cfg.histdir) )
	
        if not options.combine:
	  if os.path.isfile("%s.pevents" % tmpnam):
		q=open("%s.pevents" % tmpnam, "r")
		pevents=q.read()
	  else: 
		print "Warning, %s.pevents doesn't exist.\nEither you only intend to run on Data or something went wrong (check the log and re-run if you are running on MC)."%tmpnam
		pevents=""

	if options.runIC:
	  for i in xrange(len(files)):
		jobname =  "%s%d"%(options.outputScript,i)
		jobbasename = os.path.basename(jobname)
		f = open("%s.sh"%(jobname),"w")
		
		f.write("#!/bin/bash\n")
		f.write("source /vols/cms/grid/setup.sh\n")
		f.write("export X509_USER_PROXY=/home/hep/%s/.myproxy\n"%options.user)
		f.write("export SCRAM_ARCH=slc5_amd64_gcc462\n")
		f.write("cd %s\n"%mydir)
		f.write("touch %s.sh.run\n" % os.path.join(mydir,jobname))
		f.write("eval `scramv1 runtime -sh`\n")
                whatRun = options.combine and "combiner.py" or "fitter.py"
		whatRun += "-l %s" % options.label
		f.write("if ( python %s -i %s -n %d -j %d )\n "%(whatRun,tmpnam,int(options.nJobs),i))
		f.write("then\n")
		
		f.write("   touch %s.sh.done\n" % os.path.join(mydir,jobname))
		f.write("else\n")
		f.write("   touch %s.sh.fail\n" % os.path.join(mydir,jobname))
		f.write("fi\n")
		f.write("rm %s.sh.run\n" % os.path.join(mydir,jobname))

		
	else:
	  for i in xrange(len(files)):
		jobname =  "%s%d"%(options.outputScript,i)
		jobbasename = os.path.basename(jobname)
		f = open("%s.sh"%(jobname),"w")
		
		f.write("#!/bin/bash\n")
		f.write("set -x\n")
		
		f.write("touch %s.sh.run\n" % os.path.join(mydir,jobname))
		
		f.write("cd %s\n"%mydir)
		f.write("eval `scramv1 runtime -sh`\n")
		
		f.write("cd -\n")
		f.write("cp -p %s/../libLoopAll.so . \n" % (mydir) )
		f.write("mkdir scratch\n")
		
		f.write("tar zxfCv %s/%s.tgz scratch\n" % (mydir,options.outputScript) )
		
		f.write("cd scratch\n")
		
		f.write("cat > %s.dat << EOF\n" % jobbasename)
		f.write(datfile.replace("$input_files",files[i]).replace("$histdir",""))
		f.write("\nEOF\n")
		
                if not options.combine:
		  f.write("cat > %s.dat.pevents << EOF\n" % jobbasename)
	  	  f.write(pevents)
		  f.write("\nEOF\n")
		
		f.write("ls\n")
		whatRun = options.combine and "combiner.py" or "fitter.py"
		jobfiles = ""
		runopt = ""
		if options.watchDutyCycle:
			runopt += " --watchDutyCycle"
		if options.mountEos:
			runopt += " --mountEos"
		if i < options.nJobs:
			runopt += " -n %d -j %d" % ( int(options.nJobs),i )
			for fn in ["","histograms_"]+options.addfiles:
				jobfiles += "%s%s_%d.%s " % ( fn, cfg.histfile[0], i, cfg.histfile[1] )
			jobfiles += "%s_%d.%s " % ( cfg.outfile[0], i, cfg.outfile[1] )
			if not options.combine:
				jobfiles += "%s_%d.json " % ( cfg.histfile[0], i )
				jobfiles += "histograms_%s_%d.csv " % ( cfg.histfile[0], i )
		else:
			for fn in ["","histograms_"]+options.addfiles:
				jobfiles += "%s%s.%s:%s%s_%d.%s " % ( fn, cfg.histfile[0], cfg.histfile[1], fn, cfg.histfile[0], i, cfg.histfile[1] )
			jobfiles += "%s.%s:%s_%d.%s " % ( cfg.outfile[0], cfg.outfile[1], cfg.outfile[0], i, cfg.outfile[1] )
			if not options.combine:
				jobfiles += "%s.json:%s_%d.json " % ( cfg.histfile[0], cfg.histfile[0], i )
				jobfiles += "histograms_%s.csv:histograms_%s_%d.csv " % ( cfg.histfile[0], cfg.histfile[0], i )
		
		f.write("""
filesToCopy=\"%(files)s\"
dstFolder=\"%(histdir)s\"
cp=\"%(cp)s\"

python %(run)s -i %(job)s.dat %(runopt)s
retval=$?

if [[ $retval == 0 ]]; then
    errors=\"\"
    for f in $filesToCopy; do
         set $(echo $f | tr ':' ' ')
         $cp $1 $dstFolder/$2
         if [[ $? != 0 ]]; then
             errors=\"$errors $f($retval)\"
         fi
    done
    if [[ -z \"$errors\" ]]; then
        touch %(statfile)s.sh.done
    else
        echo 100 > %(statfile)s.sh.fail
    fi
else
    echo $retval > %(statfile)s.sh.fail
fi
"""                 % { "run":whatRun,
			"runopt":runopt,
			"files": jobfiles,
			"cp":cp,
			"job":jobbasename,
			"statfile": os.path.join(mydir,jobname),
			"histdir":cfg.histdir
			})
		f.write("rm %s.sh.run\n" % os.path.join(mydir,jobname))
		f.close()
		os.system("chmod 755 %s.sh"%(jobname))
		
	print "Submission Scripts written %sN.sh N=0,%d"%(options.outputScript,len(files)-1)

#	if not options.runIC:
	print "Written ", outnam
	print "Combine workspace after running with python combiner.py -i ",outnam

	cmdline=os.path.join(scriptdir, "cmdline.sh")
	g = open(cmdline,"w+")
	g.write("#!/bin/bash\n")
	for a in sys.argv:
		g.write("%s "%a)
	g.write("@$\n")
	g.close()
