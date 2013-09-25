import commands,sys,os

siteHandling = {
   "cern.ch"    : { "ls"      :  "xrd eoscms.cern.ch ls %s | awk '{ print $5 }' | sort",
                    "prepend" :  "root://eoscms//eos/cms" },
   "T2_CH_CSCS" : { "ls"      : "xrd cms01.lcg.cscs.ch ls %s | awk '{ print $5 }' | sort ",
                    "prepend" : "root://cms01.lcg.cscs.ch/"  }
   }


def makeCaFiles(dir,njobs=-1,jobid=0,nf=[0],maxfiles=-1,site="cern.ch"):

   dir = str(dir)
   return_files = []

   replace = None
   ls = None
   prepend = None
   if dir.startswith("/castor"):
      ls = 'nsls %s'
      prepend = 'rfio://'
   else:
      sh = siteHandling[site]
      while type(sh) == str:
         sh = siteHandling[sh]         
      ls = sh["ls"]
      prepend = sh.get("prepend",None)
      replace = sh.get("replace",None)
      
   sc,flist = commands.getstatusoutput(ls%dir)
   
   if not sc:
      files = flist.split('\n')
      ifile = 0
      for f in files:
         if '.root' in f:
            if( maxfiles > 0 and ifile >= maxfiles):
               break
            ifile += 1
            nf[0] += 1
            fname = f
            if replace:
               fname = fname.replace( *replace )
            if prepend:
               fname = "%s%s" % ( prepend, fname)
            if (njobs > 0) and (nf[0] % njobs != jobid):
               return_files.append((fname,False))
	    else:
               return_files.append((fname,True))
   else:
      sys.exit("No Such Directory: %s"%(dir))

   if nf[0]==0:
      sys.exit("No .root Files found in directory - %s:\n%s"%(dir,flist))

   return return_files


def makeDcFiles(dir,njobs=-1,jobid=0,nf=[0],maxfiles=-1):

   dcache_prepend = 'root://xrootd.grid.hep.ph.ic.ac.uk/'
   dir = str(dir)
   return_files = []

   sc,flist = commands.getstatusoutput('srmls $DCACHE_SRM_ROOT/%s --count 1000'%(dir))
   
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
   if sc==0: #status must be completed for unmount
      sc,flist = commands.getstatusoutput("rmdir %s"%dir) 
   else:
      sys.exit("Unmount "+dir+" failed. Exiting.")
   

def makeEosFiles(dir,njobs=-1,jobid=0,nf=[0]):
   sys.exit("makeEosFiles not supported anymore")
   
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

if __name__ == "__main__":
   from sys import argv
   site = "cern.ch"
   if len(argv) == 3:
      site = argv.pop(1)
   flist = makeCaFiles(argv[1],site=site)
   for f in flist:
      print f[0]
      
      
black_list = ["root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/processed/V13_03_05/data/DoublePhotonPromptReco2012B/PromptPhoton2012Data_628_1_MEr.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/processed/V13_03_05/mc/Summer12_S7_8TeV/VBF_HToGG_M-145_8TeV_sub2/SignalMC_19_2_g4J.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/processed/V13_03_06/mc/Summer12_S7_8TeV/GluGluToHToGG_M-110_8TeV/Signal_MC_3_1_QrU.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/processed/V13_03_06/mc/Summer12_S7_8TeV/GluGluToHToGG_M-120_8TeV/Signal_MC_11_1_YC2.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/HCP2012_preapproval_red_v1/data/Photon_Run2012A-recover-06Aug2012-v1__sub1/Photon_Run2012A-recover-06Aug2012-v1__sub1_0.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/hcp2012_unblind_reduction_v2/mc/Summer12_S10_8TeV/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_13.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/hcp2012_unblind_reduction_v2/mc/Summer12_S10_8TeV/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_31.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/hcp2012_unblind_reduction_v2/mc/Summer12_S10_8TeV/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_7.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_vbf_optimization_v1/mc/Summer12_RD_DR53X-PU_S10_START53_V7C/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_S10_START53_V7C-v1_v3/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_S10_START53_V7C-v1_v3_86.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/moriond2013_reduction_v1/mc/Summer12_S10_8TeV/QCD_Pt-30to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_ff/QCD_Pt-30to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_ff_23.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/VBF_HToGG_M-150_8TeV-powheg-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1/VBF_HToGG_M-150_8TeV-powheg-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1_0.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/VBF_HToGG_M-150_8TeV-powheg-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1/VBF_HToGG_M-150_8TeV-powheg-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1_7.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/VBF_HToGG_M-150_8TeV-powheg-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1/VBF_HToGG_M-150_8TeV-powheg-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1_8.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/VBF_HToGG_M-150_8TeV-powheg-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1/VBF_HToGG_M-150_8TeV-powheg-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1_9.root",
              
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_0.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_2.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_23.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_24.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_25.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_27.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_28.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_29.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_3.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_30.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_33.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_47.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_48.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_50.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_51.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_52.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_54.root",
              "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_55.root",
			"root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_56.root",
			"root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_57.root",
			"root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_58.root",
			"root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_59.root",
			"root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_60.root",
			"root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_61.root",
                       "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_62.root",
                       "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_64.root",
                       "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_65.root",
                       "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_66.root",
                       "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_67.root",
                       "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_68.root",
                       "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_69.root",
                       "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_70.root",
                       "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_71.root",
                       "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_72.root",
                       "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_73.root",
                       "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_75.root",
                       "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_76.root",
                       "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1/DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_78.root",
                    "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/reduced/reduction_RDMC_June2013_v2/mc/TTH_HToGG_M-150_8TeV-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v2/TTH_HToGG_M-150_8TeV-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v2_0.root"
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_196_1_7J4.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_316_1_wrw.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_723_1_ilG.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_804_1_Yjp.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1109_1_LoI.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1202_1_sht.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1204_1_FVO.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1225_1_vzK.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1233_1_OQk.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1236_1_dLW.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1255_1_n87.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1267_1_Vgo.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1277_1_jOY.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1283_1_o2i.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1290_1_uHE.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1291_1_TfS.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1297_1_xNF.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1308_1_sW0.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1309_1_fOn.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1311_1_NfE.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1316_1_mB2.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1324_1_U05.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1330_1_odc.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1331_1_YV2.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1336_1_yDU.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1351_1_wLJ.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1357_1_Se0.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1370_1_qom.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1378_1_sPB.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1379_1_vri.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1380_1_SJS.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1391_1_wGK.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1415_1_S4a.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1417_1_gsl.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1418_1_Me0.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1419_1_qA4.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1430_1_U7l.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1432_1_sL3.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1451_1_hfo.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1452_1_gdJ.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1477_1_Ab0.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1480_1_5vk.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1483_1_i05.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1507_1_yae.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1511_1_SaK.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1517_1_poi.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1534_1_kKk.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1543_1_xGr.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1545_1_Ymb.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1659_2_lgD.root",
                       "root://cms01.lcg.cscs.ch//store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/GJet_M80_doubleEMEnriched_8TeV-sherpa_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1813_1_jZA.root",
"root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_1_1_SxP.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_2_1_QLS.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_3_1_Oza.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_4_1_wQC.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_5_1_i2i.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_6_1_zd7.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_7_1_V81.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_8_1_cUI.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_9_1_fMt.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_10_1_WDU.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_11_1_oPO.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_12_1_KbL.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_13_1_gqT.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_14_1_E76.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_15_1_taP.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_16_1_gnQ.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_17_1_gGL.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_18_1_azV.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_19_1_YJC.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_20_1_L0u.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_21_1_U3v.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_22_1_g9M.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_23_1_uuF.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_24_1_JUf.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_25_1_vLt.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_26_1_vG1.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_27_1_zZ4.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_28_1_V0y.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_29_1_mb3.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_30_1_YqM.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_31_1_IHq.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_32_1_sqD.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_33_1_R91.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_34_1_xj5.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_35_1_wm5.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_36_1_e9S.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_37_1_Zjq.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_38_1_w9j.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_38_2_rq7.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_39_1_vvt.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_40_1_woz.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_41_1_bzt.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_42_1_Y3H.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_43_1_xuk.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_44_1_ta0.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_45_1_9ic.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_46_1_ty5.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_47_1_YV9.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_48_1_pqO.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_49_1_u79.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_50_1_Ktj.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_51_1_mS1.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_52_1_uUe.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_53_1_zU5.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_54_1_Gn5.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_55_1_ybP.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_56_1_n3E.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_57_1_NrZ.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_58_1_RzJ.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_59_1_ih4.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_59_2_A3g.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_60_1_rHZ.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_61_1_qlL.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_62_1_oU2.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_63_1_zZG.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_64_1_ODu.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_65_1_fGR.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_66_1_X7C.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_67_1_9b4.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_68_1_4Jm.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_69_1_uHY.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_70_1_ybY.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_71_1_WHX.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_72_1_xSG.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_73_1_z68.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_74_1_m8D.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_75_1_np4.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_76_1_tdS.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_77_1_u1I.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_78_1_ghg.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_79_1_p4X.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_80_1_pVh.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_81_1_xvc.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_82_1_EhC.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_83_1_FpR.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_84_1_N1B.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_85_1_iQy.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_86_1_qlF.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_87_1_gsT.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_88_1_hOe.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_89_1_bpU.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_90_1_Mjw.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_91_1_KFV.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_92_1_tkl.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_93_1_W6R.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_94_1_iSl.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_95_1_aYY.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_96_1_ZO0.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_97_1_cfw.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_98_1_3IG.root",
                       "root://cms01.lcg.cscs.ch//store/user/amarini/cmshgg/V15_00_11_v2/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM/Diphoton2Jets_EW4_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM_99_1_xGu.root",
                       "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/Summer12_DR53X-PU_RD1_START53_V7N/Graviton2PMGluGluToHToGG_M-125_8TeV-jhu-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1/Graviton2PMGluGluToHToGG_M-125_8TeV-jhu-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1_62_1_gkO.root",
                       "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/Summer12_RD1/GluGluToHToGG_M-123_8TeV-powheg-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1/GluGluToHToGG_M-123_8TeV-powheg-pythia6_Summer12_DR53X-PU_RD1_START53_V7N-v1_5_1_9v7.root",
                       "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/processed/V15_00_05/data/DoublePhoton_Run2012B-22Jan2013-v1_AOD/DoublePhoton_Run2012B-22Jan2013-v1_AOD_89_2_Is6.root",
                       ]

