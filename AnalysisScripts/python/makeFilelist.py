import commands,sys,os

def makeCaFiles(dir,njobs=-1,jobid=0,nf=[0]):

   dir = str(dir)
   return_files = []

   sc,flist = commands.getstatusoutput('nsls %s'%(dir))
#   nf = 0
   
   if not sc:
      files = flist.split('\n')
      for f in files:
         if '.root' in f:
            nf[0] += 1
            if (njobs > 0) and (nf[0] % njobs != jobid):
               return_files.append(('rfio://'+dir+'/'+f,False))
	    else:
               return_files.append(('rfio://'+dir+'/'+f,True))
            
   else:
      sys.exit("No Such Directory: %s"%(dir))

   if nf[0]==0: sys.exit("No .root Files found in directory - %s"%dir)

   return return_files


def makeDcFiles(dir,njobs=-1,jobid=0,nf=[0]):

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
         if '.root' in f:
            nf[0] += 1
            if (njobs > 0) and (nf[0] % njobs != jobid):
               return_files.append((dcache_prepend+dir+'/'+f,False))
	    else:
               return_files.append((dcache_prepend+dir+'/'+f,True))

   else:
      sys.exit("No Such Directory: %s"%(dir))
      
   if nf[0]==0: sys.exit("No .root Files found in directory - %s"%dir)

   return return_files


def makeFiles(dir,njobs=-1,jobid=0,nf=[0]):

   dir = str(dir)
   return_files = []
#   nf = 0
   if os.path.isdir(dir): 
      files = os.listdir(dir)
      for f in files:
         if '.root' in f:
            nf[0] += 1
            if (njobs > 0) and (nf[0] % njobs != jobid):
               return_files.append((dir+'/'+f,False))
	    else:
               return_files.append((dir+'/'+f,True))
   else: sys.exit("No Such Directory as %s"%dir)  

   if nf[0]==0: sys.exit("No .root Files found in directory - %s"%dir)
   
   return return_files  
