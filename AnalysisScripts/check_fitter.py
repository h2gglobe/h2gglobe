#!/usr/bin/env python

import sys, glob, os, commands
from mk_fitter import Conf 

def usage():
    print "Usage: check_fitter.py <task_folder> [post_processing_script]"
    print ""
    print " The script checks the status of all fitter jobs in <task_folder>"
    print "  Upon completion of all jobs, all workspaces are combined and the script returrns with 0."
    print "  Otherwise the script returns with -1."
    print ""
    print " A continuos monitoring of the folder can be performed with the following command line:"
    print "   while( ! ./check_fitter.py <task_folder> [post_processing_script] ); do sleep 360; done"
    print ""
    print " Optionally, a post-processing script can be called after the workspaces have been merged."
    print " The script is rung with <taks_folder> as argument."
    print "  See do_intepolation.sh for an example."
    print ""
    

if len(sys.argv) == 1:
    usage()
    sys.exit(-1)


taskdir = sys.argv[1]
postProc=None
if( len(sys.argv) > 2 ):
    postProc=sys.argv[2]

jobs = glob.glob( "%s/sub*.sh" % taskdir )
status = {}

for j in jobs:
    status[j] = ""
    for stat in "run","fail","done":
        if os.path.isfile("%s.%s" % (j,stat) ):
            status[j] = stat

groups = { "" : [], "run" : [], "fail" : [], "done" : [] } 
for i,s in status.iteritems():
    groups[s].append(i)

print taskdir
for g,jo in groups.iteritems():
    print
    print "%d jobs in status %s" % (len(jo), g)
    print "     ",
    for j in jo: print "%s" % os.path.basename(j).replace("sub","").replace(".sh",""),
    print

autorestart = [20,21,152]
restart = ""
for j in groups["fail"]:
    try:
        st=open("%s.fail" % j)
        errcode=int(st.read().split("\n")[0])
        st.close()
        if errcode in autorestart:
            jobid = os.path.basename(j).split("sub")[1].split(".")[0]
            restart += "%s," % jobid
    except:
        pass

if restart != "":
    print
    print
    print "Resubmitting jobs: %s " % restart
    os.system("./submit_fitter.py -q 8nh -j %s -d %s" % (restart, taskdir))
    print
    
if len(groups["done"]) == len(jobs):
    print "All jobs completed"
    filestocmb=glob.glob("%s/filestocombine_*.dat" % taskdir)
    print filestocmb
    if len(filestocmb) == 1:
        cfg = Conf() 
        cfile = open( filestocmb[0], "r" )

        for line in cfile.read().split("\n"):
            if line[0] == '#':
                continue
            if "histfile" in line:
                cfg.read_histfile(line)
                line = line.replace(cfg.histdir,"./")
                break
        combinedws="%s.%s" % ( os.path.join(cfg.histdir,cfg.histfile[0]), cfg.histfile[1] )
        print combinedws

        if not os.path.isfile( combinedws  ):
            os.system("python combiner.py --mountEos -i %s" % filestocmb[0] )
            if postProc:
                os.system("%s %s | tee -a %s/post_proc.log" % ( postProc, taskdir, taskdir ) )

    sys.exit(0)
            
            
sys.exit(-1)
