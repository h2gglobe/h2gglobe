{
  gSystem->Load("libPhysics.so");
  gSystem->Load("libCore.so");
  gSystem->Load("libLoopAll.so");

  gBenchmark->Start("Reduction");
  Util* ut = new Util();

  ut->SetTypeRun(1, "GluGlu2H2GG140_reduced.root");
  ut->AddFile("GluGluToHToGG_M-140_7TeV-powheg-pythia6_3_1_NqG.root",1);
  
  ut->LoopAndFillHistos();
  gBenchmark->Show("Reduction");
}

