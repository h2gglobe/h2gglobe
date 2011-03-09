{
  gSystem->Load("libPhysics.so");
  gSystem->Load("libCore.so");
  gSystem->Load("libLoopAll.so");

  gBenchmark->Start("Analysis");
  Util* ut = new Util();

  ut->SetTypeRun(2, "hist.root");
  ut->AddFile(" reducedExamples/GluGlu2H2GG140_reduced.root", 1);
  ut->AddFile(" reducedExamples/GluGlu2H2GG140_3_reduced.root", 1);
  ut->AddFile(" reducedExamples/DiPhotonBox_reduced.root", -1);
  ut->AddFile(" reducedExamples/QCD_Pt40_reduced.root", -2);
  
  ut->LoopAndFillHistos();
  gBenchmark->Show("Analysis");

  ut->WriteHist();  
}

