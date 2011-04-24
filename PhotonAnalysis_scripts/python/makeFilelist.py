
import commands,sys

def makeCaFiles(dir):

   dir = str(dir)
   return_files = []

   sc,flist = commands.getstatusoutput('nsls %s'%(dir))

   if not sc:
    files = flist.split('\n')
    for f in files:
	if '.root' in f:
	 return_files.append('rfio:'+dir+'/'+f)

   else:
    sys.exit("No Such Directory: %s"%(dir))
   return return_files

