void cloneOpttree() { 
  
  TFile *oldfile = new TFile("opttree_merged_v2.root");
  TTree *oldtree = (TTree*)oldfile->Get("opttree");
  Float_t itype, dipho_mva, cicpf4cutlevel1, cicpf4cutlevel2;

  oldtree->SetBranchAddress("itype", &itype);
  oldtree->SetBranchAddress("dipho_mva", &dipho_mva);
  oldtree->SetBranchAddress("cicpf4cutlevel1", &cicpf4cutlevel1);
  oldtree->SetBranchAddress("cicpf4cutlevel2", &cicpf4cutlevel2);
  Long64_t nentries = oldtree->GetEntries();

  TFile *newfile = new TFile("common.root", "recreate");
  TTree *newtree = oldtree->CloneTree(0);
  
  for (Long64_t z=0; z<nentries; z++) {
    oldtree->GetEntry(z);
    if (itype == 0) {
      if (dipho_mva > -0.05 && cicpf4cutlevel1>=4 && cicpf4cutlevel2>=4) {
	itype = -37;
	newtree->Fill();
      }
    }
  }

  newtree->Print();
  newtree->AutoSave();
  delete oldfile;
  delete newfile;
}
