import commands,sys,os

def makeCaFiles(dir,njobs=-1,jobid=0,nf=[0],maxfiles=-1):

   dir = str(dir)
   return_files = []

   try:
      ld_path = os.getenv("LD_LIBRARY_PATH")
   except:
      ld_path = ""

   if not "/afs/cern.ch/project/eos/installation/pro/lib64/" in ld_path:
      os.putenv("LD_LIBRARY_PATH", "%s:%s" % ( ld_path, "/afs/cern.ch/project/eos/installation/pro/lib64/" ) )

   iscastor = False
   if dir.startswith("/castor"):
      sc,flist = commands.getstatusoutput('nsls %s'%(dir))
      prepend = 'rfio://'
   else:
      ### sc,flist = commands.getstatusoutput("cmsLs %s | awk '{ print $5 }' | xargs cmsPfn | sed 's/\?.*$//'" % (dir))
      eos = "/afs/cern.ch/project/eos/installation/pro/bin/eos"
      sc,flist = commands.getstatusoutput("%s root://eoscms ls /eos/cms%s" % (eos,dir))
      prepend = 'root://eoscms//eos/cms'
      
   if not sc:
      files = flist.split('\n')
      ifile = 0
      for f in files:
         if '.root' in f:
            if( maxfiles > 0 and ifile >= maxfiles):
               break
            ifile += 1
            nf[0] += 1
            fname = "%s%s/%s" % ( prepend, dir, f)
            if (njobs > 0) and (nf[0] % njobs != jobid):
               return_files.append((fname,False))
	    else:
               return_files.append((fname,True))      
   else:
      sys.exit("No Such Directory: %s"%(dir))

   if nf[0]==0: sys.exit("No .root Files found in directory - %s"%dir)

   return return_files


def makeDcFiles(dir,njobs=-1,jobid=0,nf=[0],maxfiles=-1):

   #dcache_prepend = 'dcap://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms'
   dcache_prepend = 'root://xrootd.grid.hep.ph.ic.ac.uk/'
   dir = str(dir)
   return_files = []

   sc,flist = commands.getstatusoutput('srmls $DCACHE_SRM_ROOT/%s --count 1000'%(dir))

 #  nf = 0
   
   if not sc:
      files = flist.split('\n')
      for f in files:
         if len(f) < 1: continue
         f = (f.split()[-1]).split('/')[-1]
         ifile = 0
         if '.root' in f:
            if( maxfiles > 0 and ifile >= maxfiles):
               break
            ifile += 1
            nf[0] += 1
            if (njobs > 0) and (nf[0] % njobs != jobid):
               return_files.append((dcache_prepend+dir+'/'+f,False))
	    else:
               return_files.append((dcache_prepend+dir+'/'+f,True))

   else:
      sys.exit("No Such Directory: %s"%(dir))
      
   if nf[0]==0: sys.exit("No .root Files found in directory - %s"%dir)

   return return_files

def unmounteos(dir):
   
   unmount = 'csh -c "eosumount %s "'%dir
   sc,flist = commands.getstatusoutput(unmount) 
   sc,flist = commands.getstatusoutput("rm -r %s"%dir) 
   

def makeEosFiles(dir,njobs=-1,jobid=0,nf=[0]):
   eoslsprepend = "./eosfoo/cms/store/group/phys_higgs/cmshgg/" 
   eosprepend = "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/" 

   fc = commands.getstatusoutput('csh -c "eosmount ./eosfoo /store/cms/"') 
   dir = str(dir)
   return_files = []
#   nf = 0
   if os.path.isdir(eoslsprepend+dir): 
      files = os.listdir(eoslsprepend+dir)
      for f in files:
         ifile = 0
         if '.root' in f:
            if( maxfiles > 0 and ifile >= maxfiles):
               break
            ifile += 1
            nf[0] += 1
            if (njobs > 0) and (nf[0] % njobs != jobid):
               return_files.append((eosprepend+dir+'/'+f,False))
	    else:
               return_files.append((eosprepend+dir+'/'+f,True))
   else: 
	unmounteos("eosfoo")
	sys.exit("No Such Directory as %s"%eoslsprepend+dir)  

   if nf[0]==0: 
	unmounteos("eosfoo")
	sys.exit("No .root Files found in directory - %s"%eoslsprepend+dir)
   
   unmounteos("eosfoo")
   return return_files  

def makeFiles(dir,njobs=-1,jobid=0,nf=[0],maxfiles=-1):

   dir = str(dir)
   return_files = []
#   nf = 0
   if os.path.isdir(dir): 
      files = os.listdir(dir)
      for f in files:
         ifile = 0
         if '.root' in f:
            if( maxfiles > 0 and ifile >= maxfiles):
               break
            ifile += 1
            nf[0] += 1
            if (njobs > 0) and (nf[0] % njobs != jobid):
               return_files.append((dir+'/'+f,False))
	    else:
               return_files.append((dir+'/'+f,True))
   else: sys.exit("No Such Directory as %s"%dir)  

   if nf[0]==0: sys.exit("No .root Files found in directory - %s"%dir)
   
   return return_files  
