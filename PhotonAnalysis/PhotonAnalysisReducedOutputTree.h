void LoopAll::PhotonAnalysisReducedOutputTree() {
  utilInstance->outputTree->Branch("pho_n", &pho_n, "pho_n/I");
  utilInstance->outputTree->Branch("pho_p4", "TClonesArray", &pho_p4, 32000, 0);

  utilInstance->outputTree->Branch("pho_calopos", "TClonesArray", &pho_calopos, 32000, 0);
  utilInstance->outputTree->Branch("pho_hoe",  &pho_hoe, "pho_hoe[pho_n]/F");
  utilInstance->outputTree->Branch("pho_sieie", &pho_sieie,"pho_sieie[pho_n]/F");
  utilInstance->outputTree->Branch("pho_ecalsumetconedr03", &pho_ecalsumetconedr03, "pho_ecalsumetconedr03[pho_n]/F");
  utilInstance->outputTree->Branch("pho_ecalsumetconedr04", &pho_ecalsumetconedr04, "pho_ecalsumetconedr04[pho_n]/F");
  utilInstance->outputTree->Branch("pho_hcalsumetconedr03", &pho_hcalsumetconedr03, "pho_hcalsumetconedr03[pho_n]/F");
  utilInstance->outputTree->Branch("pho_hcalsumetconedr04", &pho_hcalsumetconedr04, "pho_hcalsumetconedr04[pho_n]/F");
  utilInstance->outputTree->Branch("pho_trksumptsolidconedr03", &pho_trksumptsolidconedr03, "pho_trksumptsolidconedr03[pho_n]/F");
  utilInstance->outputTree->Branch("pho_trksumptsolidconedr04", &pho_trksumptsolidconedr04, "pho_trksumptsolidconedr04[pho_n]/F");
  utilInstance->outputTree->Branch("pho_trksumpthollowconedr03", &pho_trksumpthollowconedr03, "pho_trksumpthollowconedr03[pho_n]/F");
  utilInstance->outputTree->Branch("pho_trksumpthollowconedr04", &pho_trksumpthollowconedr04, "pho_trksumpthollowconedr04[pho_n]/F");
  utilInstance->outputTree->Branch("pho_haspixseed", &pho_haspixseed, "pho_haspixseed[pho_n]/F");
  utilInstance->outputTree->Branch("pho_r9", &pho_r9, "pho_r9[pho_n]/F");
  utilInstance->outputTree->Branch("pho_r1x5", &pho_r1x5, "pho_r1x5[pho_n]/F");
  utilInstance->outputTree->Branch("pho_r2x5", &pho_r2x5, "pho_r2x5[pho_n]/F");
  utilInstance->outputTree->Branch("pho_scind", &pho_scind, "pho_scind[pho_n]/I");

}

