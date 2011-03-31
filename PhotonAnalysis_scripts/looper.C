{
  gSystem->Load("libPhysics.so");
  gSystem->Load("libCore.so");
  gSystem->Load("libLoopAll.so");

  gBenchmark->Start("Analysis");
  Util* ut = new Util();

  ut->ReadInput(2);
  
  ut->LoopAndFillHistos();
  gBenchmark->Show("Analysis");

  ut->WriteHist();  
  ut->WriteCounters();  
}

