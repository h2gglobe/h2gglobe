import ROOT
from python.configProducer import *

ROOT.gSystem.Load("libRooFit.so")
ROOT.gSystem.Load("libPhysics.so");
ROOT.gSystem.Load("libCore.so");
ROOT.gSystem.Load("../libLoopAll.so");

ROOT.gBenchmark.Start("Reduction")

ut = ROOT.Util();
cfg = configProducer(ut,"filestoreduce.dat",1)
ut.LoopAndFillHistos()

ROOT.gBenchmark.Show("Reduction");

#{
#  gSystem->Load("libPhysics.so");
#  gSystem->Load("libCore.so");
#  gSystem->Load("../libLoopAll.so");
#
#  gBenchmark->Start("Reduction");
#  Util* ut = new Util();
#
#  ut->ReadInput(1);
#  
#  ut->LoopAndFillHistos();
#  gBenchmark->Show("Reduction");
#}

