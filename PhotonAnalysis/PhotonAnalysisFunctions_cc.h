#define PADEBUG 0

void LoopAll::InitRealPhotonAnalysis(Util * ut, int typerun) {

  // Book histos only if not reduce step
  if (typerun != 1) {   
    for( int ind=0; ind<ut->ntypes; ind++) {
      histoContainer[ind]->Add("pho_pt", 100, 0, 100);
      histoContainer[ind]->Add("invmass_barrel", 200, 0, 200);
      histoContainer[ind]->Add("invmass_endcap", 200, 0, 200);
    }
  }
}

void LoopAll::TermRealPhotonAnalysis(int typerun) 
{}


void LoopAll::myFillHistPhotonAnalysis(Util *ut, int jentry) {

  if(PADEBUG) 
    cout << "myFillHist START"<<endl;

  counters[0]++;

  int histVal = ut->type2HistVal[ut->datatype[ut->current]];
  b_pho_n->GetEntry(jentry); 
  b_pho_p4->GetEntry(jentry); 
  
  for (int i=0; i<pho_n; i++) {
    TLorentzVector *p4 = (TLorentzVector *) pho_p4->At(i);
    histoContainer[histVal]->Fill("pho_pt", p4->Pt());
  }
  
  Int_t in_endcap = 0;
  Float_t best_mass = 0;
  for (int i=0; i<pho_n-1; i++) {
    TLorentzVector *pg1= (TLorentzVector *) pho_p4->At(i);
    if (fabs(pg1->Eta()) > 1.479)
      in_endcap = 1;

    for (int j=i+1; j<pho_n; j++) {
      TLorentzVector *pg2= (TLorentzVector *) pho_p4->At(j);
      if (fabs(pg2->Eta()) > 1.479)
	in_endcap = 1;
      TLorentzVector higgs = (*pg1) + (*pg2);
      Float_t mass = higgs.M();
      if (mass > best_mass)
	best_mass = mass;
    }
  }     

  if (best_mass != 0) 
    histoContainer[histVal]->Fill("invmass", best_mass);
  
  if(PADEBUG) 
    cout<<"myFillHist END"<<endl;
}


void LoopAll::myFillHistPhotonAnalysisRed(Util * ut, int jentry) {

  if(PADEBUG) 
    cout << "myFillHistRed START"<<endl;
  
  int histVal = ut->type2HistVal[ut->datatype[ut->current]];
  
  counters[0]++;

  b_pho_n->GetEntry(jentry); 
  b_pho_p4->GetEntry(jentry); 
  
  for (int i=0; i<pho_n; i++) {
    TLorentzVector *p4 = (TLorentzVector *) pho_p4->At(i);
    histoContainer[histVal]->Fill("pho_pt", p4->Pt());
  }
  
  Int_t in_endcap = 0;
  Float_t best_mass = 0;
  for (int i=0; i<pho_n-1; i++) {
    TLorentzVector *pg1= (TLorentzVector *) pho_p4->At(i);
    if (fabs(pg1->Eta()) > 1.479)
      in_endcap = 1;

    for (int j=i+1; j<pho_n; j++) {
      TLorentzVector *pg2= (TLorentzVector *) pho_p4->At(j);
      if (fabs(pg2->Eta()) > 1.479)
	in_endcap = 1;
      TLorentzVector higgs = (*pg1) + (*pg2);
      Float_t mass = higgs.M();
      if (mass > best_mass)
	best_mass = mass;
    }
  }     

  if (best_mass != 0) {
    if (in_endcap == 0)
      histoContainer[histVal]->Fill("invmass_barrel", best_mass);
    else
      histoContainer[histVal]->Fill("invmass_endcap", best_mass);
  }

  if(PADEBUG) 
    cout<<"myFillHistRed END"<<endl;
}

void LoopAll::myReducePhotonAnalysis(Util * ut, int jentry) {

  if(PADEBUG) 
    cout<<"myReducePhotonAnalysis START"<<endl;

  //count all events
  countersred[0]++;

  if(outputFile) {
    if(makeOutputTree) {
      
      //first selection and fill output tree
      if(!myFillReducedVarPhotonAnalysis(ut, jentry)) 
	return;
      
      //additional selection
      if(!mySelectEventRedPhotonAnalysis(ut, jentry)) 
	return;

      countersred[1]++;

      outputEvents++;
      if(PADEBUG) 
	cout<<"before fill"<<endl;

      outputTree->Fill();
      if(PADEBUG) 
	cout<<"after fill"<<endl;

      if(outputEvents==100) {
	outputEvents=0;
	outputTree->Write(0,TObject::kWriteDelete);
      }
    }
  }

  if(PADEBUG) 
    cout<<"myReducePhotonAnalysis END"<<endl;
}

void LoopAll::myGetBranchPhotonAnalysis() {
  b_pho_n = fChain->GetBranch("pho_n");
  b_pho_p4 = fChain->GetBranch("pho_p4");
}

int LoopAll::myFillReducedVarPhotonAnalysis(Util * ut, int jentry) {
  if(PADEBUG) 
    cout<<"myFillReduceVar START"<<endl;
  
   return 1;

  if(PADEBUG) 
    cout<<"myFillReduceVar END"<<endl;

}

void LoopAll::mySetBranchAddressRedPhotonAnalysis() {
  
  fChain->SetBranchAddress("pho_n", &pho_n, &b_pho_n);
  fChain->SetBranchAddress("pho_p4", &pho_p4, &b_pho_p4);
}


int LoopAll::mySelectEventRedPhotonAnalysis(Util * ut, int jentry) {
  
  // preselection at the end
  int selectevent=0;

  b_pho_n->GetEntry(jentry);

  if (pho_n > 1)
    selectevent = 1;
  else
    selectevent = 0;

  return selectevent;
}

