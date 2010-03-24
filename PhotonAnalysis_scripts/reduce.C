{
  gSystem->Load("libPhysics.so");
  gSystem->Load("libCore.so");
  gSystem->Load("libLoopAll.so");

  gBenchmark->Start("Reduction");
  Util* ut = new Util();

  ut->SetTypeRun(1, "hgg_reduced.root");
  ut->AddFile("../test/hgg.root");
  
  ut->LoopAndFillHistos();
  gBenchmark->Show("Reduction");
}

