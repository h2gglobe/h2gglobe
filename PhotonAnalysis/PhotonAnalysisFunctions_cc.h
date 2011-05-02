#define PADEBUG 0

void LoopAll::TermRealPhotonAnalysis(int typerun) 
{
//   if (typerun==3){	
//      rooContainer->FitToData("exp","mass");
//   }

}

void LoopAll::InitRealPhotonAnalysis(int typerun) {
  if(PADEBUG) 
    cout << "InitRealPhotonAnalysis START"<<endl;


  if (typerun == 2 || typerun == 1){

  }
  
  if (typerun == 3) {  
     //RooFitting type
/*
     rooContainer->AddRealVar("mass",50.,250.);
     rooContainer->AddRealVar("mu",-0.04,-1.,-0.001);

     // -------------------------------------//
     std::vector<const char*> pars(2,"t");	 
     pars[0] = "mass";
     pars[1] = "mu";
     // -------------------------------------//
     rooContainer->AddGenericPdf("exp",
	"exp((@0)*(@1))",pars,10);

     rooContainer->CreateDataSet("mass");
*/
  }
  if(PADEBUG) 
    cout << "InitRealPhotonAnalysis END"<<endl;

}

void LoopAll::myGetEntryPhotonRedAnalysis(int jentry) {
    GetEntry(inputBranches, jentry);
}

void LoopAll::myFillHistPhotonAnalysis(int jentry) {

  if(PADEBUG) 
    cout << "myFillHist START"<<endl;

  
  for (int i=0; i<pho_n; i++) {
    TLorentzVector *p4 = (TLorentzVector *) pho_p4->At(i);
    FillHist("pho_pt", p4->Pt());
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
    FillHist("invmass", best_mass);
  
  if(PADEBUG) 
    cout<<"myFillHist END"<<endl;
}


void LoopAll::myFillHistPhotonAnalysisRed(int jentry) {

  if(PADEBUG) 
    cout << "myFillHistRed START"<<endl;

  counters[0]++;

// From Here is the Standard Dec Review Selection/ gen Level studies

  std::vector<PhotonCandidate> preselected_photons;  

  TVector3 *calopos;	
  TLorentzVector *p4;
  for (int i=0; i<pho_n; i++) {
    p4 = (TLorentzVector *) pho_p4->At(i);
    calopos  = (TVector3 *) pho_calopos->At(i);
    float pt  = p4->Pt(); 
    float eta = fabs(calopos->Eta());
    //PreSelection
     if ( 
       (! pho_haspixseed[i])
       && pt > 30. 
       && pho_hoe[i] <  0.1
       && pho_trksumptsolidconedr03[i] < 2*(3.5 + 0.001*pt)
       && pho_ecalsumetconedr03[i] < 2*(4.2 + 0.006*pt)
       && pho_hcalsumetconedr03[i] < 2*(2.2 + 0.0025*pt)
       &&((eta < 1.4442) || ((eta > 1.566) && (eta < 2.5))) 
       ) {
         PhotonCandidate candidate;
         candidate.p4 		= p4;
	 candidate.calopos	= calopos;
         candidate.pixSeed 	= pho_haspixseed[i];
         candidate.trkIso 	= pho_trksumpthollowconedr04[i];
         candidate.ecalIso 	= pho_ecalsumetconedr04[i];
         candidate.hcalIso 	= pho_hcalsumetconedr04[i];
         candidate.sieie 	= pho_sieie[i];
         candidate.hoe 		= pho_hoe[i];
         candidate.r9 		= pho_r9[i];
         preselected_photons.push_back(candidate);
       }
  }

  //Event Selection
  int n_preselected_pho = preselected_photons.size();
  //FillHist("h_n_sel",n_preselected_pho);

  // Sort Photons into Pt Order
  std::sort(preselected_photons.begin()
           ,preselected_photons.end()
           ,PhoP4greater); 
 
// Regular Event Selection begins here
  float best_mass = 0.;
  float best_pt = -1;
  int category = -1;
  float min_r9;
  float max_eta;

  if (n_preselected_pho > 1 ){

    PhotonCandidate leading, nleading;


    leading  = preselected_photons[0];
    nleading = preselected_photons[1];


  	if (leading.p4->Pt() > 40.){

      // Determine the Category of the event
      // -> Which histogram is filled
      min_r9  = min(leading.r9, nleading.r9);
	    max_eta = max(fabs(leading.calopos->Eta())
			              ,fabs(nleading.calopos->Eta()));
	 
      if (min_r9 < 0.93 && max_eta < 1.4442 ) category = 1;
      if (min_r9 > 0.93 && max_eta < 1.4442 ) category = 2;
      if (min_r9 < 0.93 && max_eta > 1.566 && max_eta < 2.5) category = 3;
      if (min_r9 > 0.93 && max_eta > 1.566 && max_eta < 2.5) category = 4;

      // -------------------------------------------------------
      TLorentzVector Higgs = (*(preselected_photons[0].p4))
                                 +(*(preselected_photons[1].p4));
           
      float mass = Higgs.M();
      float h_pt = Higgs.Pt();
      if (mass > 100. && mass < 150.){
        //Good event, passes preselection and acceptance cuts

        int pass_selection[2];
        int pass_isolation[2];
        int in_iso_gap[2];
    
        //Now do selection on leading photon
        pass_selection[0] = leading.hoe < 0.02
           	       && (((leading.sieie < 0.01)  && (fabs(leading.calopos->Eta()) < 1.4442)) 
           		  || (( leading.sieie < 0.028)
              	       && ((fabs(leading.calopos->Eta()) < 2.5) && (fabs(leading.calopos->Eta()) > 1.566))) );
        pass_isolation[0] =  leading.trkIso < (1.5 + 0.001*leading.p4->Pt())		
           	       && leading.ecalIso < (2.0 + 0.006*leading.p4->Pt())
                          && leading.hcalIso < (2.0 + 0.0025*leading.p4->Pt());
	 
        //Selection on next to leading photon
        pass_selection[1] = nleading.hoe < 0.02
           	       && (((nleading.sieie < 0.01)  && (fabs(nleading.calopos->Eta()) < 1.4442)) 
           		  || (( nleading.sieie < 0.028)
             	       && ((fabs(nleading.calopos->Eta()) < 2.5) && (fabs(nleading.calopos->Eta()) > 1.566))) );
        pass_isolation[1] =  nleading.trkIso < (1.5 + 0.001*nleading.p4->Pt())
           	       && (nleading.ecalIso < (2.0 + 0.006*nleading.p4->Pt()))
           	       && (nleading.hcalIso < (2.0 + 0.0025*nleading.p4->Pt()));

              
	        if (pass_selection[0] && pass_isolation[0] ){

	            if (pass_selection[1] && pass_isolation[1]){
        	      cout << "mass is " << mass << " and higgs pt is " << h_pt << endl;
		          FillHist("pho_pt",category,leading.p4->Pt());
		          FillHist("pho_pt",category,nleading.p4->Pt());
              		  best_mass = mass;
 		          best_pt   = h_pt;
            }
	      }
      }
    }
  }

  if (best_mass > 100)
    FillHist("mass",0, best_mass);
  FillHist("pt",0, best_pt);
  if (category > -1){
    FillHist("mass",category, best_mass);
    FillHist("pt",category, best_pt);
  }
  if(PADEBUG) 
    cout<<"myFillHistRed END"<<endl;
}

void LoopAll::myStatPhotonAnalysis(int jentry) {

  if(PADEBUG) 
    cout << "myStat START"<<endl;
/*
  counters[0]++;
     
  std::vector<PhotonCandidate> preselected_photons;  

  TVector3 *calopos;	
  TLorentzVector *p4;
  for (int i=0; i<pho_n; i++) {
    p4 = (TLorentzVector *) pho_p4->At(i);
    calopos  = (TVector3 *) pho_calopos->At(i);
    float pt  = p4->Pt(); 
    float eta = fabs(calopos->Eta());
    //PreSelection
     if ( 
       (! pho_haspixseed[i])
       && pt > 30. 
       && pho_hoe[i] <  0.1
       && pho_trksumptsolidconedr03[i] < 2*(3.5 + 0.001*pt)
       && pho_ecalsumetconedr03[i] < 2*(4.2 + 0.006*pt)
       && pho_hcalsumetconedr03[i] < 2*(2.2 + 0.0025*pt)
       &&((eta < 1.4442) || ((eta > 1.566) && (eta < 2.5))) 
       ) {
         PhotonCandidate candidate;
         candidate.p4 		= p4;
	 candidate.calopos	= calopos;
         candidate.pixSeed 	= pho_haspixseed[i];
         candidate.trkIso 	= pho_trksumpthollowconedr04[i];
         candidate.ecalIso 	= pho_ecalsumetconedr04[i];
         candidate.hcalIso 	= pho_hcalsumetconedr04[i];
         candidate.sieie 	= pho_sieie[i];
         candidate.hoe 		= pho_hoe[i];
         preselected_photons.push_back(candidate);
       }
  }

  //Event Selection
  int n_preselected_pho = preselected_photons.size();

  if (n_preselected_pho > 1 ){

     // Sort Photons into Pt Order
     std::sort(preselected_photons.begin()
           ,preselected_photons.end()
           ,PhoP4greater); 
 
    // Regular Event Selection begins here

     PhotonCandidate leading, nleading;

     leading  = preselected_photons[0];
     nleading = preselected_photons[1];


  	 if (leading.p4->Pt() > 40.){
         

         // -------------------------------------------------------
         TLorentzVector Higgs = (*(preselected_photons[0].p4))
                                 +(*(preselected_photons[1].p4));
           float mass = Higgs.M();
             if (mass > 50. && mass < 200.){
             //Good event, passes preselection and acceptance cuts

             int pass_selection[2];
             int pass_isolation[2];
             int in_iso_gap[2];
    
            //Now do selection on photons
             if(		 leading.hoe < 0.02
                	       && (((leading.sieie < 0.01)  && (fabs(leading.calopos->Eta()) < 1.4442)) 
                		  || (( leading.sieie < 0.028)
                   	       && ((fabs(leading.calopos->Eta()) < 2.5) && (fabs(leading.calopos->Eta()) > 1.566))) )
             		       && leading.trkIso < (1.5 + 0.001*leading.p4->Pt())		
                	       && leading.ecalIso < (2.0 + 0.006*leading.p4->Pt())
                               && leading.hcalIso < (2.0 + 0.0025*leading.p4->Pt())
	 			
             		       && nleading.hoe < 0.02
                	       && (((nleading.sieie < 0.01)  && (fabs(nleading.calopos->Eta()) < 1.4442)) 
                		  || (( nleading.sieie < 0.028)
                  	       && ((fabs(nleading.calopos->Eta()) < 2.5) && (fabs(nleading.calopos->Eta()) > 1.566))) )
             		       && nleading.trkIso < (1.5 + 0.001*nleading.p4->Pt())
                	       && (nleading.ecalIso < (2.0 + 0.006*nleading.p4->Pt()))
                	       && (nleading.hcalIso < (2.0 + 0.0025*nleading.p4->Pt()))
	        )
		rooContainer->SetRealVar("mass",mass);
               }
	}
    }
  
*/
}


void LoopAll::myReducePhotonAnalysis(int jentry) {

  if(PADEBUG) 
    cout<<"myReducePhotonAnalysis START"<<endl;

  //count all events
  countersred[0]++;

  if(outputFile) {
    if(makeOutputTree) {
      
      //first selection and fill output tree
      if(!myFillReducedVarPhotonAnalysis(jentry)) 
	return;
      
      //additional selection
      if(!mySelectEventRedPhotonAnalysis(jentry)) 
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
  /// GetBranches(inputBranchNames, inputBranches);
}


void LoopAll::mySetBranchAddressRedPhotonAnalysis() {
  std::cerr << "mySetBranchAddressRedPhotonAnalysis" << std::endl;
  inputBranches.clear();
  GetBranches(inputBranchNames, inputBranches);
  SetBranchAddresses(inputBranchNames);
}


int LoopAll::myFillReducedVarPhotonAnalysis(int jentry) {
  if(PADEBUG) 
    cout<<"myFillReduceVar START"<<endl;
  
   return 1;

  if(PADEBUG) 
    cout<<"myFillReduceVar END"<<endl;

}

int LoopAll::mySelectEventRedPhotonAnalysis(int jentry) {
  
  // preselection at the end
  int selectevent=0;

  b_pho_n->GetEntry(jentry);

  if (pho_n > 1)
    selectevent = 1;
  else
    selectevent = 0;

  return selectevent;
}

