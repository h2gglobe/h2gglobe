void LoopAll::PhotonAnalysisReducedOutputTree() {
  UtilInstance->outputTree->Branch("pho_n", &pho_n, "pho_n/I");
  UtilInstance->outputTree->Branch("pho_p4", "TClonesArray", &pho_p4, 32000, 0);
}

