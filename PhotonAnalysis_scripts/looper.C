{
  gSystem->Load("libPhysics.so");
  gSystem->Load("libCore.so");
  gSystem->Load("libLoopAll.so");

  gBenchmark->Start("Analysis");
  Util* ut = new Util();

  ut->SetTypeRun(2, "hist.root");
  ut->AddFile("hgg_reduced.root");
  
  ut->LoopAndFillHistos();
  gBenchmark->Show("Anslysis");

  ut->WriteHist();  
}

