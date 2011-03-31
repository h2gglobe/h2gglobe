{
  gSystem->Load("libPhysics.so");
  gSystem->Load("libCore.so");
  gSystem->Load("libLoopAll.so");

  gBenchmark->Start("Reduction");
  Util* ut = new Util();

  ut->ReadInput(1);
  
  ut->LoopAndFillHistos();
  gBenchmark->Show("Reduction");
}

