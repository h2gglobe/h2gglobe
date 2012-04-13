#!/usr/bin/env python
import sys, os, commands

class Conf:
	def __init__(self):
		pass

	def __str__(self):
		rt = ""
		for k,v in self.__dict__.iteritems():
			rt += "%s = %s\n" % (k,v)
		return rt
	
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
	parser.add_option("-i","--inputDat",dest="inputDat")
	parser.add_option("-n","--nJobs",type="int", dest="nJobs",default=-1)
	parser.add_option("-o","--outputScript",dest="outputScript",default="sub")
	parser.add_option("-l","--label",dest="label",default="")
	(options,args)=parser.parse_args()
	
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

		if in_comment:
			datfile += "%s\n" % line
			continue
		
		line = line % { "label" : options.label } 

		if "histfile" in line:
			cfg.read_histfile(line)
			line = line.replace(cfg.histdir,"./")
			

		if line.startswith("split"):
			files.append("%s\n" % line.replace("split ",""))
			if not has_polder:
				datfile += "%s"
				has_polder = True
		elif "typ" in line:
			for i in xrange(options.nJobs):
				files[i] += "%s\n" % line
			if not has_polder:
				datfile += "%s"
				has_polder = True
		else:
			datfile += "%s\n" % line

	mydir=os.getcwd()
	scriptdir=os.path.dirname(options.outputScript)
	os.system("mkdir -p %s" % scriptdir)

	outnam=os.path.join(scriptdir, "filestocombine_%s" % os.path.basename(options.inputDat))
	g = open(outnam,"w+") 
	filestocmb = ""
	for i in xrange(len(files)):
		fil = commands.getoutput("cmsPfn %s_%d.%s" % ( os.path.join(cfg.histdir,cfg.histfile[0]), i, cfg.histfile[1] ))
		filestocmb += "Fil=%s\n" % fil
	g.write( (datfile % filestocmb).replace("histfile=./","histfile=%s" % scriptdir ) )
	g.close()	

	tmpnam = os.path.join(os.path.dirname(options.inputDat), "tmp_%s" % os.path.basename(options.inputDat))
	if not os.path.isfile("%s.pevents" % tmpnam):
		print "Generating the pevents file...",
		tmp  = open("%s" % tmpnam, "w+")
		idat = open(options.inputDat, "r")
		tmp.write(idat.read().replace("split",""))
		tmp.close()
		idat.close()
		print "python fitter.py -i %s --dryRun >& %s.log\n" % (tmpnam,tmpnam)
		os.system("python fitter.py -i %s --dryRun >& %s.log\n" % (tmpnam,tmpnam) )
		print "Done. Check %s.log for errors" % tmpnam
		
	mkdir="mkdir"
	cp="cp -pv"
	if cfg.histdir.startswith("/castor"):
		mkdir="rfmkdir"
		cp="rfcp"
	elif cfg.histdir.startswith("/store"):
		cp="cmsStage"
		mkdir="cmsMkdir"
		

	os.system("tar zcf %s.tgz $(find -name \*.dat -or -name \*.py) aux" % options.outputScript)
	os.system("%s %s" % ( mkdir, cfg.histdir) )
	os.system("%s %s.tgz %s" % ( cp,  options.outputScript, cfg.histdir) )
	
	q=open("%s.pevents" % tmpnam, "r")
	pevents=q.read()

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
		f.write(datfile % files[i])
		f.write("\nEOF\n")
		
		f.write("cat > %s.dat.pevents << EOF\n" % jobbasename)
		f.write(pevents)
		f.write("\nEOF\n")
		
		f.write("ls\n")
		
		if i < options.nJobs:
			f.write("if ( python fitter.py -i %s.dat -n %d -j %d ) "%(jobbasename,int(options.nJobs),i))
			for fn in "","histograms_":
				f.write("&& ( %s %s%s_%d.%s %s ) "        % ( cp, fn, cfg.histfile[0], i, cfg.histfile[1], cfg.histdir ) )
			#f.write("&& ( %s %s_%d.%s_ascii_events.txt %s ) " % ( cp, cfg.histfile[0], i, cfg.histfile[1], cfg.histdir ) )
			f.write("&& ( %s %s_%d.%s %s ) " % ( cp, cfg.outfile[0], i, cfg.outfile[1], cfg.histdir ) )
			f.write("&& ( %s %s_%d.json %s ) "                % ( cp, cfg.histfile[0], i, cfg.histdir ) )
			f.write("&& ( %s histograms_%s_%d.csv %s ) " % ( cp, cfg.histfile[0], i, cfg.histdir ) )
		else:
			f.write("if ( python fitter.py -i %s.dat ) "%(jobbasename))
			for fn in "","histograms_":
				f.write("&& ( %s %s%s.%s %s/%s%s_%d.%s ) " % ( cp, fn, cfg.histfile[0], cfg.histfile[1], cfg.histdir, fn, cfg.histfile[0], i, cfg.histfile[1] ) )
			#f.write("&& ( %s %s.%s_ascii_events.txt %s/%s_%d.%s_ascii_events.txt ) " % ( cp, cfg.histfile[0], cfg.histfile[1], cfg.histdir, cfg.histfile[0], i, cfg.histfile[1] ) )
			f.write("&& ( %s %s.%s %s/%s_%d.%s ) " % ( cp, cfg.outfile[0], cfg.outfile[1], cfg.histdir, cfg.outfile[0], i, cfg.outfile[1]) )
			f.write("&& ( %s %s.json %s/%s_%d.json ) " % ( cp, cfg.histfile[0], cfg.histdir, cfg.histfile[0], i ) )
			f.write("&& ( %s histograms_%s.csv %s/%s_%d.csv ) " % ( cp, cfg.histfile[0], cfg.histdir, cfg.histfile[0], i ) )
			
				
		f.write("; then\n")
		
		f.write("   touch %s.sh.done\n" % os.path.join(mydir,jobname))
		f.write("else\n")
		f.write("   touch %s.sh.fail\n" % os.path.join(mydir,jobname))
		f.write("fi\n")
		
		f.write("rm %s.sh.run\n" % os.path.join(mydir,jobname))
		f.close()
		os.system("chmod 755 %s.sh"%(jobname))
		
	print "Submission Scripts written %sN.sh N=0,%d"%(options.outputScript,len(files)-1)
	print "Written ", outnam
	print "Combine workspace after running with python combiner.py -i ",outnam
