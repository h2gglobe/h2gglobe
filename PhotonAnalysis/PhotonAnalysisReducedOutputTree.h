void LoopAll::PhotonAnalysisReducedOutputTree() {
  utilInstance->outputTree->Branch("pho_n", &pho_n, "pho_n/I");
  utilInstance->outputTree->Branch("pho_p4", "TClonesArray", &pho_p4, 32000, 0);
}

