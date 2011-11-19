//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//==============================  Interactive Plotting  ===============================
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#define MPDEBUG 1

void LoopAll::myPlotInteractive(TString hsmallname) {

  if(MPDEBUG) cout<<"myPlotInteractive Start real "<<endl;
  cout<<"myPlotInteractive Start real "<<endl;

  TFile *file1;

  if(hsmallname!="") {
    file1 = TFile::Open(hsmallname);
    file1->cd();
  }

  if(MPDEBUG) cout<<"myPlotInteractive 1 "<<endl;
  cout<<"myPlotInteractive 1 "<<endl;

  //MARCOMARCO
  lineWidth = 4;
  lineWidth = 2;
  NReBin=2;
  BigLegendInd=1;
  StartCat=0;

  TH1F * hist1D[NIND];
  TH2F * hist2D[NIND];
  TProfile * histProf[NIND];
  TH1F * histSgnl;
  TH1F * histSvN;
  TH1F * histBkgd;
  TH1F * histBkgdFit;
  TH1F * histSoverB;
  TH2F * histSvsSoB;
  TH2F * histSvsB;
  TLegend *legend;
  TString histnamebase;
  TString histnamebaseSig;
  TString histnamebaseData;

  if(MPDEBUG) cout<<"myPlotInteractive 2 "<<endl;

  // OPTIONS string for prompting user
  OPTIONS = "\n\nOptions:\n";
  OPTIONS = "navigating:\n";
  OPTIONS += "\tls \t\t(list hists)\n";
  OPTIONS += "\tq \t\t(quit)\n";
  OPTIONS += "\t# \t\t(advance # histograms)\n";
  OPTIONS += "\tn \t\t(next)\n";
  OPTIONS += "\tp \t\t(preceeding)\n";
  OPTIONS += "\tb \t\t(back to last hist)\n";
  OPTIONS += "\tj #\t\t(jump to hist #)\n\n";
  OPTIONS += "manipulating histograms:\n";
  OPTIONS += "\tscale \t\t(normalize to 1 on/off)\n";
  OPTIONS += "\tstack \t\t(on/off)\n";
  OPTIONS += "\tpopsig \t\t(plot signal stack in front of bkgd)\n";
  OPTIONS += "\tomit \t\t(on/off for each file)\n";
  OPTIONS += "\trebin N\t\t(rebin in N bins)\n";
  OPTIONS += "\tint \t\t(plot integral)\n";
  OPTIONS += "\trint \t\t(plot reverse integral)\n";
  OPTIONS += "\tsb \t\t(sgnl,bkgd,and s/b on new canvi)\n";
  OPTIONS += "\tsvn\t\tplot S/sqrt(S+B)\n";
  OPTIONS += "\tsumw2 \t\t(Sumw2 on/off)\n";
  OPTIONS += "\tcats \t\t(plot all categories on/off)\n";
  OPTIONS += "\tmulti N M\t(plot next N*M hists on 1 canvas)\n";
  OPTIONS += "\tfill \t\t(on/off)\n";
  OPTIONS += "\tfit xmin xmax \t\t(fit S/B - not working yet)\n\n";
  OPTIONS += "manupulating canvas:\n";
  OPTIONS += "\tresize x y\t(set canvas size)\n";
  OPTIONS += "\tlog \t\t(on/off)\n";
  OPTIONS += "\toflow \t\t(over-flow)\n";
  OPTIONS += "\tuflow \t\t(under-flow)\n";
  OPTIONS += "\tuof \t\t(both-flows)\n";
  OPTIONS += "\tgrid \t\t(X&Y on/off)\n";
  OPTIONS += "\tgridx \t\t(on/off)\n";
  OPTIONS += "\tgridy \t\t(on/off)\n\n";
  OPTIONS += "\tlegend\t\t(on/off)\n";
  OPTIONS += "\tbigleg\t\t(make 1 legend big)\n";
  OPTIONS += "\tsetmax\t\t(set max in hist, 0 gives default)\n";
  OPTIONS += "\tsetmin\t\t(set min in hist, useful for log)\n";
  OPTIONS += "\trmtitles\t(removes the titles from all plots)\n";
  OPTIONS += "\ttr/tl/br/bl\t(move legend to top-right/top-left/...)\n\n";
  OPTIONS += "other:\n";
  OPTIONS += "\tpf  \t\t(print under/over-flows)\n";
  OPTIONS += "\tpb  \t\t(print bin values)\n";
  OPTIONS += "\twrite (prompts for extension to write hist to file)\n";
  OPTIONS += "\twriteall (prompts for extension to write all hists to directory allfiles/)\n";
  OPTIONS += "\twritesob (prompts for extension to write s/b hist to file)\n";

  std::pair<int,int> ipair(0,0);
  //********************************************

  ch = new TCanvas("ch","ch",10,10,xSize,ySize);
  int ivar=ipair.first;
  int icat=ipair.second;
  nrows=1;
  ncolumns=1;
  if(MPDEBUG) {
    for( int iInd=0;iInd < NIND;++iInd) {
      //std::cout << "mp->infoind[iInd]" << mp->infoind[iInd] << std::endl;
      //std::cout << "mp->itypePI[mp->infoind[iInd]]" << mp->itypePI[mp->infoind[iInd]] << std::endl;
    }
  }
  if(MPDEBUG) cout<<"myPlotInteractive 3 "<<endl;


  


  do {
    ivar=ipair.first;
    icat=ipair.second;


    if(MPDEBUG) cout<<"myPlotInteractive 4 ivar icat: " << ivar << " " << icat << endl;
    if(DoRedraw) {
      Int_t varDim=1;
      //if(h2d[ivar]==1)varDim=2;
      varDim=h2d[ivar]+1;
      if(typplot[ivar]==-1) typplot[ivar]=typplotall;
      //cout <<"Plotting histo n. "<<ivar<<" "<<h2d[ivar]<<"  "<<plotvarnames[ivar]<<endl;
      int ncat=histoncat[ivar];
      if(ncat==0) {
        ncat=1;
      }
      if(ncat==1) {
        nrows=1;
        ncolumns=1;
      }
      if(MPDEBUG)std::cout << "plotvarnames[ivar]: " << plotvarnames[ivar] << std::endl;
      if(ncat>1 && DoCats) {
        if (ncat==2)    {  nrows=1;      ncolumns=2; }
        else if (ncat==3) {  nrows=1;      ncolumns=3; }
        else if (ncat<5)  {  nrows=2;      ncolumns=2; }
        else if (ncat<7)  {  nrows=2;      ncolumns=3; }
        else if (ncat<9)  {  nrows=2;      ncolumns=4; }
        else if (ncat<10) {  nrows=3;      ncolumns=3; }
        else if (ncat<13) {  nrows=4;      ncolumns=3; }
        else if (ncat<17) {  nrows=4;      ncolumns=4; }
        else              {  nrows=5;      ncolumns=5; }
      }
      if(MPDEBUG) cout<<"myPlotInteractive 5 "<<endl;
      
      //   if(!DoReSize){
      //   xSize = singlexSize*nrows;
      //   ySize = singleySize*ncolumns;
      //   if(xSize>xSizeMAX)xSize=xSizeMAX;
      //   if(ySize>ySizeMAX)ySize=ySizeMAX;
      //   }
       
      ch->Clear();
      ch->SetWindowSize(xSize,ySize);
      ch->SetCanvasSize(xSize-4,ySize-28);
      if(DoSnB) {
        ch2->Clear();
        ch2->SetWindowSize(xSize,ySize);
        ch2->SetCanvasSize(xSize-4,ySize-28);
        ch3->Clear();
        ch3->SetWindowSize(xSize,ySize);
        ch3->SetCanvasSize(xSize-4,ySize-28);
      }
      if(DoSvB) {
        ch4->Clear();
        ch4->SetWindowSize(xSize,ySize);
        ch4->SetCanvasSize(xSize-4,ySize-28);
        ch5->Clear();
        ch5->SetWindowSize(xSize,ySize);
        ch5->SetCanvasSize(xSize-4,ySize-28);
      }
      if(MPDEBUG) cout<<"myPlotInteractive 6 "<<endl;
      Int_t EndCat = ncat;
      if(StartCat>=ncat) {
        StartCat=0;
      } else if(!DoCats) {
        EndCat=TMath::Min(ncat,StartCat+nrows*ncolumns);
      }
      if(DoBigLegend && (BigLegendInd < StartCat || BigLegendInd > EndCat))BigLegendInd=1;




      for(Int_t jcat = StartCat; jcat != EndCat;++jcat) {

        //set canvas attributes
        if(ncat>1 && DoMulti) {
          if(jcat==StartCat)ch->Divide(ncolumns,nrows);
          if(DoSnB) {
            if(jcat==StartCat)ch2->Divide(ncolumns,nrows);
            if(jcat==StartCat)ch3->Divide(ncolumns,nrows);
            ch2->cd(jcat-StartCat+1);
            ch2->cd(jcat-StartCat+1)->SetLogy(DoLog);
            ch2->cd(jcat-StartCat+1)->SetGridx(DoGridX);
            ch2->cd(jcat-StartCat+1)->SetGridy(DoGridY);
            ch3->cd(jcat-StartCat+1);
            ch3->cd(jcat-StartCat+1)->SetLogy(DoLog);
            ch3->cd(jcat-StartCat+1)->SetGridx(DoGridX);
            ch3->cd(jcat-StartCat+1)->SetGridy(DoGridY);
          }
          if(DoSvB) {
            if(jcat==0)ch4->Divide(ncolumns,nrows);
            ch4->cd(jcat-StartCat+1);
            //ch4->cd(jcat-StartCat+1)->SetLogy(DoLog);
            ch4->cd(jcat-StartCat+1)->SetGridx(DoGridX);
            ch4->cd(jcat-StartCat+1)->SetGridy(DoGridY);
            if(jcat==0)ch5->Divide(ncolumns,nrows);
            ch5->cd(jcat-StartCat+1);
            //ch5->cd(jcat-StartCat+1)->SetLogy(DoLog);
            ch5->cd(jcat-StartCat+1)->SetGridx(DoGridX);
            ch5->cd(jcat-StartCat+1)->SetGridy(DoGridY);
          }
          ch->cd(jcat-StartCat+1);
          ch->cd(jcat-StartCat+1)->SetLogy(DoLog);
          ch->cd(jcat-StartCat+1)->SetGridx(DoGridX);
          ch->cd(jcat-StartCat+1)->SetGridy(DoGridY);
        } else { 
          if(DoSnB) {
            ch2->cd();
            ch2->SetLogy(DoLog);
            ch2->SetGridx(DoGridX);
            ch2->SetGridy(DoGridY);
            ch3->cd();
            ch3->SetLogy(DoLog);
            ch3->SetGridx(DoGridX);
            ch3->SetGridy(DoGridY);
          }
          if(DoSvB) {
            ch4->cd();
            //ch4->SetLogy(DoLog);
            ch4->SetGridx(DoGridX);
            ch4->SetGridy(DoGridY);
            ch5->cd();
            //ch5->SetLogy(DoLog);
            ch5->SetGridx(DoGridX);
            ch5->SetGridy(DoGridY);
          }
          ch->cd();
          ch->SetLogy(DoLog);
          ch->SetGridx(DoGridX);
          ch->SetGridy(DoGridY);
        }
        if(!DoMulti && jcat != icat)continue;

        if(MPDEBUG) cout<<"myPlotInteractive 7 "<<endl;
        //place legend
        Float_t BLF = 1.;// Big Legend Factor
        if(DoBigLegend && jcat==BigLegendInd-1)BLF=1.5;
        Float_t IBLF = 1./BLF;// Inverted Big Legend Factor
        if(legendx1>0&&legendy1>0&&legendx2>legendx1&&legendy2>legendy1) {
          legend = new TLegend(legendx1,legendy1,legendx2,legendy2);
        } else {
          if(DoLegendT  &&  DoLegendR)legend = new TLegend(IBLF*0.7,IBLF*0.7,0.89,0.89);
          if(DoLegendT  && !DoLegendR)legend = new TLegend(0.11,IBLF*0.7,BLF*0.3,0.89);
          if(!DoLegendT &&  DoLegendR)legend = new TLegend(IBLF*0.7,0.11,0.89,BLF*0.3);
          if(!DoLegendT && !DoLegendR)legend = new TLegend(0.11,0.11,BLF*0.3,BLF*0.3);
        }
        if(MPDEBUG) cout<<"myPlotInteractive 71 "<<endl;
        legend->SetBorderSize(0);
        legend->SetTextFont(62);
        legend->SetLineColor(0);
	//MARCOMARCO
        legend->SetLineStyle(1);
        legend->SetLineWidth(lineWidth);
        legend->SetFillColor(kWhite);
        legend->SetFillStyle(0);
        if(MPDEBUG) cout<<"myPlotInteractive 72 "<<endl;

        //get hist names
        histnamebase=plotvarnames[ivar];

        if(MPDEBUG) cout<<"myPlotInteractive 73 "<<endl;
        TString catstring="_cat"+num_to_string(jcat);
        if(xaxislabels[ivar] == "xaxis") {
          if(histoncatindtonames[ivar]>=0)catstring="_"+plotvarcatnames[histoncatindtonames[ivar]][jcat];
          if(ncat>0)histnamebase+=catstring;
        } else {
          if(histoncatindtonames[ivar]>=0)catstring=plotvarcatnames[histoncatindtonames[ivar]][jcat];
          if(ncat>0)histnamebase+=catstring;
        }
        if(MPDEBUG) cout<<"myPlotInteractive 74 "<<endl;
        histnamebaseSig=histnamebase+"Sig";
        histnamebaseData=histnamebase+"Data";
        if(MPDEBUG) cout<<"myPlotInteractive 75 "<<endl;


        if(MPDEBUG) cout<<"myPlotInteractive 75 "<<NIND<<endl;

        for(Int_t iInd=0;iInd!=NIND;++iInd){
        if(MPDEBUG) cout<<"myPlotInteractive 76 "<<omitInd[iInd]<<endl;
          if(omitInd[iInd]%2==1)continue;
          ////myHistname[iInd]=plotvarnames[ivar]+"_"+num_to_string(mp->histoind[iInd]);
          //myHistname[iInd]=plotvarnames[ivar]+"_"+num_to_string(histoind[infoind[iInd]]);
          myHistname[iInd]=plotvarnames[ivar]+"_"+inshortnames[iInd];
          if(ncat>0) {
            ////myHistname[iInd]=plotvarnames[ivar]+"_cat"+num_to_string(jcat)+"_"+num_to_string(mp->histoind[iInd]);
            //myHistname[iInd]=plotvarnames[ivar]+"_cat"+num_to_string(jcat)+"_"+num_to_string(histoind[infoind[iInd]]);
            myHistname[iInd]=plotvarnames[ivar]+"_cat"+num_to_string(jcat)+"_"+inshortnames[iInd];
          }


          if(MPDEBUG)std::cout << "myHistname["<<iInd<<"]: " << myHistname[iInd] << std::endl;

	  std::cout<<"NEWNEW "<<iInd<<" "<<infoind[iInd]<<" "<<inshortnames[iInd]<<std::endl;





        }

        if(MPDEBUG) cout<<"myPlotInteractive 8 "<<endl;
        if(varDim==1) {
          if(MPDEBUG) cout<<"myPlotInteractive 9 "<<endl;
          Float_t MAXVAL = 0; 
          //Get Histograms
          float upscalemax = 1.1;
          if(DoLog)upscalemax=3.;
          for(Int_t iInd=0;iInd!=NIND;++iInd){
            if(omitInd[iInd]%2==1)continue;
            //hist1D[iInd] = new TH1F();
            //*(hist1D[iInd]) = *((TH1F*)gROOT->FindObject(myHistname[iInd])); 

            if(MPDEBUG) std::cout << "myHistname[iInd]: " << myHistname[iInd] << std::endl;

            hist1D[iInd] = new TH1F(*((TH1F*)gROOT->FindObject(myHistname[iInd]))); 
            if(MPDEBUG) std::cout << "1 myHistname[iInd]: " << myHistname[iInd] << std::endl;
            hist1D[iInd]->SetLineStyle(0);
            hist1D[iInd]->GetXaxis()->SetLabelSize(0.04);
            hist1D[iInd]->GetYaxis()->SetLabelSize(0.04);
            if(MPDEBUG) std::cout << "2 myHistname[iInd]: " << myHistname[iInd] << std::endl;
            if(xaxislabels[ivar] != "xaxis") {
              hist1D[iInd]->GetXaxis()->SetTitle(xaxislabels[ivar]);
              //std::cout <<"xaxislabels[ivar]: " << xaxislabels[ivar] << std::endl;
            }
            if(MPDEBUG) std::cout << "3 myHistname[iInd]: " << myHistname[iInd] << std::endl;
            hist1D[iInd]->GetXaxis()->SetTitleSize(0.06);
            hist1D[iInd]->GetXaxis()->SetTitleOffset(0.75);
            if(MPDEBUG) std::cout << "4 myHistname[iInd]: " << myHistname[iInd] << std::endl;
          }
          if(MPDEBUG) cout<<"myPlotInteractive 10 "<<endl;
          //reset hist contents (rebin, integrate, ...)
          for(Int_t iInd=0;iInd!=NIND;++iInd){
            if(omitInd[iInd]%2==1)continue;
            if(DoUnderFlow){
              float uflow = hist1D[iInd]->GetBinContent(0);
              float bin1content = hist1D[iInd]->GetBinContent(1);
              hist1D[iInd]->SetBinContent(1,bin1content+uflow);
            }
            if(DoOverFlow){
              Int_t binN = hist1D[iInd]->GetNbinsX();
              float oflow = hist1D[iInd]->GetBinContent(binN+1);
              float binNcontent = hist1D[iInd]->GetBinContent(binN);
              hist1D[iInd]->SetBinContent(binN,binNcontent+oflow);
            }
            if(DoRebin) hist1D[iInd]->Rebin(NReBin); 
          }
          if(MPDEBUG) cout<<"myPlotInteractive 11 "<<endl;
          if(DoSnB || DoSvB || DoSvN) {
            //signal and background
            histSgnl = new TH1F();
            histBkgd = new TH1F();
            histBkgdFit = new TH1F();
            histSoverB = new TH1F();
            histSvsB = new TH2F();
            histSvsSoB = new TH2F();
            histSvN = new TH1F();
            bool firstsgnl=true;
            bool firstbkgd=true;
            for(Int_t iInd=0;iInd!=NIND;++iInd){
              if(omitInd[iInd]%2==1)continue;
              if(mySgnl[iInd]) {
                if(firstsgnl) {
                  *(histSgnl) = *(hist1D[iInd]);
                } else {
                  histSgnl->Add(hist1D[iInd]);
                }
                firstsgnl=false;
              } else if(myBkgd[iInd]) {
                if(firstbkgd) {
                  *(histBkgd) = *(hist1D[iInd]);
                } else {
                  histBkgd->Add(hist1D[iInd]);
                }
                firstbkgd=false;
              }
            }
            if(xaxislabels[ivar] != "xaxis") {
              histBkgd->GetXaxis()->SetTitle(xaxislabels[ivar]);
              histSgnl->GetXaxis()->SetTitle(xaxislabels[ivar]);
            }
            histBkgd->SetLineColor(kBlue);
            histSgnl->SetLineColor(kRed);
          }
          // now have signal and background hists (histSgnl and histBkgd)
          if(DoSnB) {
            if(DoMulti && ncat>1) {
              ch2->cd(jcat-StartCat+1);
            } else {
              ch2->cd();
            }
            *(histBkgdFit) = *(histBkgd);
            TF1 * fsig = new TF1("fsig","gaus");
            TF1 * fback = new TF1("fback","[0]+[1]*pow(x,1)+[2]*pow(x,2)");
            if(DoFitSig) histSgnl->Fit("fsig","fitsig","",fitXmin,fitXmax);
            if(DoFit || DoFitSig) {
              histBkgdFit->Fit(fback,"","",histBkgd->GetBinCenter(1),histBkgd->GetBinCenter(histBkgd->GetNbinsX()));
              for(Int_t ibin=1;ibin<=histBkgdFit->GetNbinsX();++ibin) {
                histBkgdFit->SetBinContent(ibin,fback->Eval(histBkgdFit->GetBinCenter(ibin)));
              }
            }

            if(!DoLog) {
              histBkgd->SetMinimum(0);
	      cout<<"1 histBkgd->Draw();"<<endl;
              histBkgd->Draw();
              if(RmTitles) histBkgd->SetTitle("");
	      cout<<"2 histSgnl->Draw(same);"<<endl;
              histSgnl->Draw("same");
            } else {
              histSgnl->SetMaximum(upscalemax*histBkgd->GetBinContent(histBkgd->GetMaximumBin()));
              if(sethistmax>0)histSgnl->SetMaximum(sethistmax);
              if(sethistmin>0)histSgnl->SetMinimum(sethistmin);
	      cout<<"3 histSgnl->Draw();"<<endl;
              histSgnl->Draw();
              if(RmTitles) histBkgd->SetTitle("");
              histBkgd->Draw("same");
	      cout<<"4 histBkgd->Draw(same);"<<endl;
            }
            if(DoFit || DoFitSig) fback->Draw("same");
            std::cout << "Background at M=100GeV: " << fback->Eval(100.) << std::endl;
            std::cout << "Background at M=110GeV: " << fback->Eval(110.) << std::endl;
            std::cout << "Background at M=120GeV: " << fback->Eval(120.) << std::endl;
            std::cout << "Background at M=130GeV: " << fback->Eval(130.) << std::endl;
            if(DoFitSig)fsig->Draw("same");

            if(DoMulti && ncat>1) {
              ch3->cd(jcat-StartCat+1);
            } else {
              ch3->cd();
            }
            *(histSoverB) = *(histSgnl);

            if(DoFit) {
              histSoverB->Divide(histBkgdFit);
            } else {
              histSoverB->Divide(histBkgd);
            }
            if(!DoLog)histSoverB->SetMinimum(0);
            if(!DoCats || jcat==StartCat) {
              sobmax = histSoverB->GetBinContent(histSoverB->GetMaximumBin());
            }
            histSoverB->SetMaximum(upscalemax*sobmax);
            if(sethistmax>0)histSoverB->SetMaximum(sethistmax);
            if(sethistmin>0)histSoverB->SetMinimum(sethistmin);
            histSoverB->Draw("hist");

            if(DoFitSig) {
              histSoverB->Fit("gaus","","",fitXmin,fitXmax);
            }
            std::cout << "cat " << jcat << "\tmax S/B\t" << sobmax << std::endl;
          }
          if(DoSvB) {
            // signal vs background
            TH1F * hS = new TH1F();
            TH1F * hB = new TH1F();
            *(hS) = *(histSgnl);
            *(hB) = *(histBkgd);
            // if DoUnderFlow/DoOverFlow are true this is already done above.
            if(!DoUnderFlow) {
              float uflow = hS->GetBinContent(0);
              float bin1content = hS->GetBinContent(1);
              hS->SetBinContent(1,bin1content+uflow);
              uflow = hB->GetBinContent(0);
              bin1content = hB->GetBinContent(1);
              hB->SetBinContent(1,bin1content+uflow);
            }
            if(!DoOverFlow) {
              Int_t binN = hS->GetNbinsX();
              float oflow = hS->GetBinContent(binN+1);
              float binNcontent = hS->GetBinContent(binN);
              hS->SetBinContent(binN,binNcontent+oflow);
              binN = hB->GetNbinsX();
              oflow = hB->GetBinContent(binN+1);
              binNcontent = hB->GetBinContent(binN);
              hB->SetBinContent(binN,binNcontent+oflow);
            }
            hS->SetContent(hS->GetIntegral());
            hB->SetContent(hB->GetIntegral());
            bool flipsvb = hS->GetMean() > hB->GetMean();
            if(flipsvb) {
              for(Int_t i=0; i!=hS->GetNbinsX(); ++i) { 
                hS->SetBinContent(i+1,1.-hS->GetBinContent(i+1));
                hB->SetBinContent(i+1,1.-hB->GetBinContent(i+1));
              }
            }
            hS->GetMinimum();
            Float_t sbmin = TMath::Min(hS->GetMinimum(),hB->GetMinimum());
            Float_t sobmin = 100.;
            Float_t sobmax = 0.;
            Float_t xbkgd[hB->GetNbinsX()-2];
            Float_t ysig[hS->GetNbinsX()-2];
            Float_t xsob[hS->GetNbinsX()-2];
            for(Int_t ibin=2;ibin<hS->GetNbinsX();++ibin) {
              xbkgd[ibin-2]=hB->GetBinContent(ibin);
              xsob[ibin-2]=hS->GetBinContent(ibin)/hB->GetBinContent(ibin);
              ysig[ibin-2]=hS->GetBinContent(ibin);
              if(xsob[ibin-2]>sobmax)sobmax=xsob[ibin-2];
              if(xsob[ibin-2]<sobmin)sobmin=xsob[ibin-2];
            }
            TH2F * histSvsBtemp = new TH2F("histSvBt","Signal vs Background",hS->GetNbinsX(),0.9*sbmin,1,hS->GetNbinsX(),0.9*sbmin,1);
            TH2F * histSvsSoBtemp = new TH2F("histSvSoBt","S eff vs (S eff) / (B eff)",hS->GetNbinsX(),0.0*sobmin,1.1*sobmax,hS->GetNbinsX(),0.9*sbmin,1);
            *(histSvsB) = *(histSvsBtemp);
            *(histSvsSoB) = *(histSvsSoBtemp);
            for(Int_t ibin=1;ibin<=hS->GetNbinsX();++ibin) {
              histSvsB->Fill(hB->GetBinContent(ibin),hS->GetBinContent(ibin));
              histSvsSoB->Fill(hS->GetBinContent(ibin)/hB->GetBinContent(ibin),hS->GetBinContent(ibin));
            }
            if(DoMulti && ncat>1) {
              ch4->cd(jcat-StartCat+1);
            } else {
              ch4->cd();
            }
            TGraph * gsvsob = new TGraph(hS->GetNbinsX()-2,xsob,ysig);
            histSvsSoB->SetMarkerStyle(kStar);
            histSvsSoB->Draw();
            gsvsob->SetLineWidth(2);
            gsvsob->Draw();
            if(DoMulti && ncat>1) {
              ch5->cd(jcat-StartCat+1);
            } else {
              ch5->cd();
            }
            TGraph * gsvb = new TGraph(hS->GetNbinsX()-2,xbkgd,ysig);
            histSvsB->SetMarkerStyle(kStar);
            //histSvsB->SetMarkerColor(kWhite);
            histSvsB->Draw();
            gsvb->SetLineWidth(2);
            gsvb->Draw();
            TF1 * fline = new TF1("fline","x");
            fline->SetLineColor(kRed);
            fline->SetLineStyle(kDashed);
            fline->SetLineWidth(0);
            fline->Draw("same");
          }
          if(DoMulti && ncat>1) {
            ch->cd(jcat-StartCat+1);
          } else {
            ch->cd();
          }
          if(MPDEBUG) cout<<"myPlotInteractive 12 "<<endl;
          // scale hist
          for(Int_t iInd=0;iInd!=NIND;++iInd){
            if(omitInd[iInd]%2==1)continue;
            if(DoInt) { hist1D[iInd]->SetContent(hist1D[iInd]->GetIntegral()); }
            if(DoRevInt) {
              Int_t nx = hist1D[iInd]->GetNbinsX();
              for(Int_t i=0; i<=nx; ++i) { hist1D[iInd]->SetBinContent(i+1,1.-hist1D[iInd]->GetBinContent(i+1));}
              //Double_t *in = hist1D[iInd]->GetIntegral();
              //for(Int_t i=0; i<=nx; ++i) { hist1D[iInd]->SetBinContent(i+1,in[nx-i]); }
            }
            if(!DoInt && !DoRevInt) {
              if(!DoScale && hist1D[iInd]->Integral()>0.) {
                hist1D[iInd]->Scale(1./hist1D[iInd]->Integral());
              } else {
                if(sigScale!=0) {
                  //if(mp->itypePI[mp->infoind[iInd]]>0)hist1D[iInd]->Scale(sigScale);
                  //FIXME if(itypePI[infoind[iInd]]>0)hist1D[iInd]->Scale(sigScale);
                  if(itypePI[iInd]<0)hist1D[iInd]->Scale(sigScale);
                } else {
                  hist1D[iInd]->Scale(myScale[iInd]);
                }
              }
            } 
            Float_t tempMAXVAL = hist1D[iInd]->GetMaximum();
            if(tempMAXVAL > MAXVAL) {
              MAXVAL = tempMAXVAL;
            }
          }
          if(MPDEBUG) cout<<"myPlotInteractive 13 "<<endl;
          //set max and don't zero suppress
          if(DoInt || DoRevInt) {
            for(Int_t iInd=0;iInd!=NIND;++iInd) { if(omitInd[iInd]%2==1)continue;hist1D[iInd]->SetMaximum(1.05);}
          } else {
          float tempmax = upscalemax*MAXVAL;
          float tempmin = 0;
            if(sethistmax>0)tempmax = sethistmax;
            if(sethistmin>0)tempmin = sethistmin;
            for(Int_t iInd=0;iInd!=NIND;++iInd) { if(omitInd[iInd]%2==1)continue;hist1D[iInd]->SetMaximum(tempmax);}
            for(Int_t iInd=0;iInd!=NIND;++iInd) { if(omitInd[iInd]%2==1)continue;hist1D[iInd]->SetMinimum(tempmin);}
          }
          if(!DoLog)for(Int_t iInd=0;iInd!=NIND;++iInd) {if(omitInd[iInd]%2==1)continue; hist1D[iInd]->SetMinimum(0);}

          if(MPDEBUG) cout<<"myPlotInteractive 14 "<<endl;

          //color and fill-color hist
          for(Int_t iInd=0;iInd!=NIND;++iInd){
            if(omitInd[iInd]%2==1)continue;
            hist1D[iInd]->SetLineWidth(0); 
            if(DoFill) {
              hist1D[iInd]->SetFillColor(myColor[iInd]);
              if(omitInd[iInd]==0) {
                hist1D[iInd]->SetLineColor(kBlack); 
              } else {
                hist1D[iInd]->SetLineColor(myColor[iInd]);
              }
              //FIXME if(DoData && itypePI[infoind[iInd]]==0) {
              if(DoData && itypePI[iInd]==0) {
                hist1D[iInd]->SetFillStyle(0); 
                hist1D[iInd]->SetLineColor(kBlack); 
                //gStyle->SetErrorX(0);

		//cout<<"MARCO here marker stile and size"<<endl;
		//myMarker[iInd]=20;//marco

                hist1D[iInd]->SetMarkerStyle(myMarker[iInd]); 
                //Chris 
                //hist1D[iInd]->SetMarkerSize(1.4); 
		//if(myMarker[iInd]==20) hist1D[iInd]->SetMarkerSize(2); 

                //if(hist1D[iInd]->GetNbinsX()>50)hist1D[iInd]->SetMarkerSize(1.); 
                //if(hist1D[iInd]->GetNbinsX()>120)hist1D[iInd]->SetMarkerSize(0.5); 
              }
              //if(!DoFillSig && mp->itypePI[mp->infoind[iInd]]>0)
              //FIXME if(!DoFillSig && itypePI[infoind[iInd]]>0) {
	      if(!DoFillSig && itypePI[iInd]<0) {
                hist1D[iInd]->SetFillStyle(0); 
                hist1D[iInd]->SetLineColor(myColor[iInd]); 
                if(false && omitInd[iInd]==0)
                  hist1D[iInd]->SetLineWidth(0); 
                else 
                  hist1D[iInd]->SetLineWidth(0); 
              }
            } else {
              hist1D[iInd]->SetFillStyle(0); 
              hist1D[iInd]->SetLineColor(myColor[iInd]); 
            }
          }

          //make the stack
          if(DoStack) { 
            NStackLines=0;
            myStack = new THStack(histnamebase,histnamebase);
            mySigStack = new THStack(histnamebaseSig,histnamebaseSig);
            myDataStack = new THStack(histnamebaseData,histnamebaseData);
            if(DoPopSig) {
            bool firstbg = true;
            bool firstsig=true;
              if(MPDEBUG) cout<<"myPlotInteractive 14 - 01"<<endl;
              bool newline=true;
              Int_t nbinstemp;
              Float_t lowedgetemp;
              Float_t highedgetemp;
              for(Int_t iInd=0;iInd!=NIND;++iInd){

		cout<<"MARCO FIXME0 "<<iInd<<" "<<myStackOrder[iInd]<<endl;

                if(iInd==0) {
                  nbinstemp = hist1D[myStackOrder[iInd]]->GetNbinsX();
                  lowedgetemp = hist1D[myStackOrder[iInd]]->GetBinLowEdge(1);
                  highedgetemp = hist1D[myStackOrder[iInd]]->GetBinLowEdge(nbinstemp) + hist1D[myStackOrder[iInd]]->GetBinWidth(nbinstemp);

		  cout<<"MARCO FIXME1 "<<nbinstemp<<" "<<lowedgetemp<<" "<<highedgetemp<<endl;
                }
                if(omitInd[myStackOrder[iInd]]%2==1)continue; 
		cout<<"MARCO FIXME2 "<<nbinstemp<<" "<<lowedgetemp<<" "<<highedgetemp<<endl;

                //std::cout << "myStackOrder["<<iInd<<"]: " << myStackOrder[iInd] << "\tmp->itypePI[myStackOrder[iInd]]: " << mp->itypePI[myStackOrder[iInd]] << "\tmp->infoind[myStackOrder[iInd]]: " << mp->infoind[myStackOrder[iInd]] << std::endl;

		//MARCO FIXME
		cout<<"MARCO FIXME3 "<<iInd<<" "<<myStackOrder[iInd]<<" "<<infoind[myStackOrder[iInd]]<<endl;

		cout<<"MARCO FIXME4 "<<DoData <<" "<<itypePI[myStackOrder[iInd]]<<endl;

                if(DoData && itypePI[myStackOrder[iInd]]==0) {
                  if(MPDEBUG) cout<<"myPlotInteractive 14 - 02"<<endl;
                  std::cout << iInd << std::endl;
                  std::cout << myLabel[iInd] << std::endl;
                  std::cout << myStackOrder[iInd] << std::endl;
                  std::cout << omitInd[myStackOrder[iInd]] << std::endl;
		  std::cout << hist1D[myStackOrder[iInd]] << std::endl;
                  myDataStack->Add(hist1D[myStackOrder[iInd]]); 



		  hist1D[myStackOrder[iInd]]->Draw();

		  cout<<"drawn"<<endl;





                } else if(itypePI[myStackOrder[iInd]]<0) {
                  if(newline) {
                    if(firstsig || NStackLines==0) {
                      histStackLine[NStackLines] = TH1F("histline"+NStackLines,"",nbinstemp,lowedgetemp,highedgetemp);
                      firstsig=false;
                    } else {
                      histStackLine[NStackLines] = histStackLine[NStackLines-1];
                    }
                    histStackLineIsSig[NStackLines]=1;
                    histStackLine[NStackLines].SetFillColor(myColor[myStackOrder[iInd]]);
                    newline=false;
                  }
                  if(MPDEBUG) cout<<"myPlotInteractive 14 - 03"<<endl;
		  //MARCOMARCO
                  //mySigStack->Add(hist1D[myStackOrder[iInd]]); 
                  histStackLine[NStackLines].Add(hist1D[myStackOrder[iInd]]);
                  histStackLineColor[NStackLines] = myColor[iInd];
                  if(omitInd[myStackOrder[iInd]]==0) {
                    NStackLines++;
                    newline=true;
                  }
                } else if(itypePI[myStackOrder[iInd]]>0) {
                  if(newline) {
                    if(firstbg || NStackLines==0) {
                      histStackLine[NStackLines] = TH1F("histline"+NStackLines,"",nbinstemp,lowedgetemp,highedgetemp);
                      firstbg=false;
                    } else {
                      histStackLine[NStackLines] = histStackLine[NStackLines-1];
                    }
                    histStackLineIsSig[NStackLines]=0;
                    newline=false;
                  }
                  if(MPDEBUG) cout<<"myPlotInteractive 14 - 04"<<endl;
                  myStack->Add(hist1D[myStackOrder[iInd]]); 
                  //std::cout << "NStackLines: " << NStackLines << std::endl;
                  //std::cout << "entries: " << histStackLine[NStackLines].GetEntries() << std::endl;
                  //std::cout << nbinstemp << "\t" << lowedgetemp << "\t" << highedgetemp << std::endl;
                  histStackLine[NStackLines].Add(hist1D[myStackOrder[iInd]]);
                  histStackLineColor[NStackLines] = myColor[iInd];
                  //std::cout << "entries AFTER: " << histStackLine[NStackLines].GetEntries() << std::endl;
                  if(omitInd[myStackOrder[iInd]]==0) {
                    NStackLines++;
                    newline=true;
                  }
              if(MPDEBUG) cout<<"myPlotInteractive 14 - 04 - 01"<<endl;
                }
              if(MPDEBUG) cout<<"myPlotInteractive 14 - 04 - 02"<<endl;
              }
              if(MPDEBUG) cout<<"myPlotInteractive 14 - 04 - 03"<<endl;
            } else {  
              if(MPDEBUG) cout<<"myPlotInteractive 14 - 05"<<endl;
              bool newline=true;
              Int_t nbinstemp;
              Float_t lowedgetemp;
              Float_t highedgetemp;
              for(Int_t iInd=0;iInd!=NIND;++iInd){
                if(iInd==0) {
                  nbinstemp = hist1D[myStackOrder[iInd]]->GetNbinsX();
                  lowedgetemp = hist1D[myStackOrder[iInd]]->GetBinLowEdge(1);
                  highedgetemp = hist1D[myStackOrder[iInd]]->GetBinLowEdge(nbinstemp) + hist1D[myStackOrder[iInd]]->GetBinWidth(nbinstemp);
                }
                if(omitInd[myStackOrder[iInd]]%2==1)continue; 
                //std::cout << "POP SIG BETTER BE FALSE!!" << std::endl;
                //std::cout << "mp->itypePI[myStackOrder[iInd]]: " << mp->itypePI[myStackOrder[iInd]] << std::endl;
                if(DoData && itypePI[myStackOrder[iInd]]==0) {
                  if(MPDEBUG) cout<<"myPlotInteractive 14 - 06"<<endl;
                  myDataStack->Add(hist1D[myStackOrder[iInd]]); 
                } else {
                  if(newline) {
                    if(NStackLines==0) {
                      histStackLine[NStackLines] = TH1F("histline"+NStackLines,"",nbinstemp,lowedgetemp,highedgetemp);
                    } else {
                      histStackLine[NStackLines] = histStackLine[NStackLines-1];
                    }
                    histStackLineIsSig[NStackLines]=0;
                    newline=false;
                  }
                  if(MPDEBUG) cout<<"myPlotInteractive 14 - 07"<<endl;
                  myStack->Add(hist1D[myStackOrder[iInd]]); 
                  histStackLine[NStackLines].Add(hist1D[myStackOrder[iInd]]);
                  if(omitInd[myStackOrder[iInd]]==0) {
                    NStackLines++;
                    newline=true;
                  }
                }
              }
            }// END if(DoPopSig)
            if(MPDEBUG) cout<<"myPlotInteractive 14 - 08"<<endl;
          }// END if(DoStack)

          if(MPDEBUG) cout<<"myPlotInteractive 15 "<<endl;

          cout << "Plotting Histogram: " << ivar<<":"<<jcat << "\t" << histnamebase << endl;
          outname = "hist_"; outname += plotvarnames[ivar];//histnamebase;
          outnamesob = outname+"SoverB";

          if(MPDEBUG) cout<<"myPlotInteractive 16 "<<endl;


          //draw hist
          TString histlike = " hist ";
          TString histlikesame = " hist same ";
          if(DoSumw2)histlike="";
          if(DoSumw2)histlikesame="same";
          if(DoStack) {
            double stackmax=1;
            double stackmin=0;
            if(DoPopSig) {
              stackmax = TMath::Max(myStack->GetMaximum(),mySigStack->GetMaximum());
              if(sethistmax>0)stackmax=sethistmax;
              if(sethistmin>0)stackmin=sethistmin;
              myStack->SetMaximum(stackmax);
              mySigStack->SetMaximum(stackmax);
              myStack->SetMinimum(stackmin);
              mySigStack->SetMinimum(stackmin);
            }
            if(DoData) {
              stackmax = TMath::Max(myDataStack->GetMaximum(),stackmax);
              if(sethistmax>0)stackmax=sethistmax;
              if(sethistmin>0)stackmin=sethistmin;
              myStack->SetMaximum(stackmax);
              mySigStack->SetMaximum(stackmax);
              myDataStack->SetMaximum(stackmax);
              myStack->SetMinimum(stackmin);
              mySigStack->SetMinimum(stackmin);
              myDataStack->SetMinimum(stackmin);
            }
            if(DoData) {
              stackmax = TMath::Max(myDataStack->GetMaximum(),stackmax);
              if(sethistmax>0)stackmax=sethistmax;
              if(sethistmin>0)stackmin=sethistmin;
              myStack->SetMaximum(stackmax);
              mySigStack->SetMaximum(stackmax);
              myDataStack->SetMaximum(stackmax);
              myStack->SetMinimum(stackmin);
              mySigStack->SetMinimum(stackmin);
              myDataStack->SetMinimum(stackmin);
            }
            bool stackfilled=false;
            TList * histlist = myStack->GetHists();
            for(int ih=0;ih!=histlist->GetEntries();++ih) {
              if(((TH1*)histlist->At(ih))->Integral()>0.)stackfilled=true;
            }
            cout << "stackfilled = " << stackfilled << endl;
            //if(stackfilled) {
            cout << "99 drawing here 1" << endl; 
            myStack->Draw(histlike);
            if(xaxislabels[ivar] != "xaxis")myStack->GetXaxis()->SetTitle(xaxislabels[ivar]);
            if(yaxislabels[ivar] != "yaxis")myStack->GetYaxis()->SetTitle(yaxislabels[ivar]);
            if(RmTitles) myStack->SetTitle("");
            cout << "991 drawing here 1" << endl; 
            myStack->Draw(histlike);
            cout << "992 drawing here 1" << endl; 
            if(DoPopSig)mySigStack->Draw(histlikesame);
            //} else {
            //  if(DoPopSig) {
            //    cout << "drawing here 2" << endl; 
            //    if(MPDEBUG) cout<<"myPlotInteractive 16.51"<<endl;
            //    //if(RmTitles) myStack->SetTitle("");
            //    //myStack->Draw(histlike);
            //    if(MPDEBUG) cout<<"myPlotInteractive 16.52"<<endl;
            //    mySigStack->Draw(histlike);
            //    if(MPDEBUG) cout<<"myPlotInteractive 16.53"<<endl;
            //    //if(xaxislabels[ivar] != "xaxis")mySigStack->GetXaxis()->SetTitle(xaxislabels[ivar]);
            //    if(MPDEBUG) cout<<"myPlotInteractive 16.54"<<endl;
            //    //if(yaxislabels[ivar] != "yaxis")mySigStack->GetYaxis()->SetTitle(yaxislabels[ivar]);
            //    if(MPDEBUG) cout<<"myPlotInteractive 16.55"<<endl;
            //    //if(RmTitles) myStack->SetTitle("");
            //    if(MPDEBUG) cout<<"myPlotInteractive 16.56"<<endl;
            //    mySigStack->Draw(histlike);
            //  } else {
            //    cout << "drawing here 3" << endl; 
            //    myStack->Draw(histlike);
            //    if(xaxislabels[ivar] != "xaxis")myStack->GetXaxis()->SetTitle(xaxislabels[ivar]);
            //    if(yaxislabels[ivar] != "yaxis")myStack->GetYaxis()->SetTitle(yaxislabels[ivar]);
            //    if(RmTitles) myStack->SetTitle("");
            //    myStack->Draw(histlike);
            //  }
            //}
            if(DoData) {
              if(myStack->GetHists() || mySigStack->GetHists()){
                if(RmTitles) myDataStack->SetTitle("");
		cout << "993 drawing here 1" << endl; 
                myDataStack->Draw("same");
              } else {
                if(RmTitles) myDataStack->SetTitle("");
		cout << "994 drawing here 1" << endl; 
                myDataStack->Draw();
                if(xaxislabels[ivar] != "xaxis")myDataStack->GetXaxis()->SetTitle(xaxislabels[ivar]);
                if(yaxislabels[ivar] != "yaxis")myDataStack->GetYaxis()->SetTitle(yaxislabels[ivar]);
                if(RmTitles) myDataStack->SetTitle("");
		cout << "995 drawing here 1" << endl; 
                myDataStack->Draw();
              }
            }
//            std::cout << "NStackLines: " << NStackLines << std::endl;
//            for(int iline=0;iline!=NStackLines;++iline) {
//              std::cout << "iline: " << iline << std::endl;
//              //std:cout << "nbins: " << histStackLine[iline].GetNbinsX() << std::endl;
//              histStackLine[iline].SetLineColor(kBlack);
//              histStackLine[iline].SetLineWidth(2);
//              histStackLine[iline].Draw("same");
//            }
          } else if(DoSvN) {
            *(histSvN) = *(histSgnl);
            if(DoInt) {
              histSgnl->SetContent(histSgnl->GetIntegral());
              histBkgd->SetContent(histBkgd->GetIntegral());
            }
            if(DoRevInt) {
              Double_t *inS = histSgnl->GetIntegral();
              Double_t *inB = histBkgd->GetIntegral();
              Int_t nx = histSgnl->GetNbinsX();
              for(Int_t i=0; i<=nx; ++i) { 
                histSgnl->SetBinContent(i+1,1.-histSgnl->GetBinContent(i+1));
                histBkgd->SetBinContent(i+1,1.-histBkgd->GetBinContent(i+1));
                //histSgnl->SetBinContent(i+1,inS[nx-i]); 
                //histBkgd->SetBinContent(i+1,inB[nx-i]); 
              }
            }
            for(Int_t i=0; i!=histSgnl->GetNbinsX(); ++i) { 
              if( histBkgd->GetBinContent(i+1)+histSgnl->GetBinContent(i+1) == 0 )
                histSvN->SetBinContent(i+1,0);
              else 
                histSvN->SetBinContent(i+1,histSgnl->GetBinContent(i+1)/sqrt(histBkgd->GetBinContent(i+1)+histSgnl->GetBinContent(i+1)) );
            }
            cout << "996 drawing here 1" << endl; 
            histSvN->Draw(histlike);
          } else {
            bool firstind = true;


	    cout<<"MMM firstind NIND "<<firstind<<" "<<NIND<<endl;
            for(Int_t iInd=0;iInd!=NIND;++iInd){
              if(omitInd[iInd]%2==1)continue;
              hist1D[iInd]->SetTitle(histnamebase);
              hist1D[iInd]->SetStats(DoStats);
              if(firstind) {
                //FIXME if(DoData && itypePI[infoind[iInd]]==0) {
                if(DoData && itypePI[iInd]==0) {
                  //hist1D[iInd]->Draw("E1");
                } else {
		  cout << "997 drawing here 1" << endl; 
                  hist1D[iInd]->Draw(histlike);
                  firstind=false;
                }
              } else {
                //FIXME if(DoData && itypePI[infoind[iInd]]==0) {
                if(DoData && itypePI[iInd]==0) {
                  //hist1D[iInd]->Draw("E1 same");
                } else {
		  cout << "998 drawing here 1" << endl; 
                  hist1D[iInd]->Draw(histlike+"same");
                }
              }
            }
	    cout<<"MMM firstind NIND "<<firstind<<" "<<NIND<<endl;
            for(Int_t iInd=0;iInd!=NIND;++iInd){
              if(omitInd[iInd]%2==1)continue;
              //FIXME if(DoData && itypePI[infoind[iInd]]==0) {
              if(DoData && itypePI[iInd]==0) {
                if(firstind) {
		  cout << "999 drawing here 1" << endl; 
                  hist1D[iInd]->Draw("E1");
		  firstind=false;
                } else {
		  cout << "9910 drawing here 1" << endl; 
                  hist1D[iInd]->Draw("E1 same");
                }
              }
            }// end for(iInd)
          }

          if(MPDEBUG) cout<<"myPlotInteractive 17 "<<endl;
          //add hist to legend
          for(Int_t iInd=0;iInd!=NIND;++iInd){
            Int_t reverseInd = myStackOrder[NIND-1-iInd];// makes top of legend = top of stack
            if(omitInd[reverseInd]%2==1)continue;
            if(omitInd[reverseInd]==2)continue;
            std::cout << "myLabel: " << myLabel[reverseInd] << std::endl;

	    //marcomarco
	    if(iInd!=NIND-1) {
            if(DoFill)  legend->AddEntry(hist1D[reverseInd],myLabeldisplay[reverseInd],"f"); 
            else      legend->AddEntry(hist1D[reverseInd],myLabeldisplay[reverseInd],"l");
	    }
	    else 
	      legend->AddEntry(hist1D[reverseInd],myLabeldisplay[reverseInd],"ep");
          }
          if(MPDEBUG) cout<<"myPlotInteractive 18 "<<endl;
#ifndef MARCODEFINE
          //if(NIND > 5)legend->SetNColumns(2);// can use for newer versions of root
          if(NIND > 5)legend->SetNColumns(NColsLegend);// can use for newer versions of root
#endif
          if(DoPrintFlows || DoPrintBins) {
            for(Int_t iInd=0;iInd!=NIND;++iInd){
              if(omitInd[iInd]%2==1)continue;
              Int_t NBins = hist1D[iInd]->GetNbinsX();
              float underflow = hist1D[iInd]->GetBinContent(0);
              float overflow = hist1D[iInd]->GetBinContent(NBins+1);
              float integral=hist1D[iInd]->Integral();
              if(!DoUnderFlow)integral+=underflow;
              if(!DoOverFlow)integral+=overflow;
              if(underflow>0 || overflow>0) {
                std::cout << "\n" << myLabel[iInd] << std::endl;
                float uefficiency=0;
                float oefficiency=0;
                if(integral>0) {
                  uefficiency=underflow/integral;
                  oefficiency=overflow/integral;
                }
                if(underflow>0)std::cout << "\tunderflow\t" << setw(9) << underflow << "\t("<< uefficiency <<")"<<std::endl;
                if(overflow>0)std::cout  << "\toverflow\t"  << setw(9) << overflow  << "\t("<< oefficiency  <<")"<<std::endl;
              }
              if(DoPrintBins) {
                if(underflow==0 && overflow==0)std::cout << "\n" << myLabel[iInd] << std::endl;
                for(Int_t iBin=1;iBin<NBins+1;++iBin) {
                  float bincontent=hist1D[iInd]->GetBinContent(iBin);
                  float binefficiency=0;
                  if(integral>0)binefficiency=bincontent/integral;
                  std::cout << myLabel[iInd] << " bin: " << iBin << "\t" << setw(9) << bincontent << "\t("<<binefficiency<<")"<<std::endl;
                }
              }
              //std::cout << myLabel[iInd] << "\ttotal\t" << hist1D[iInd]->Integral() << std::endl;
            }
          }//end if(DoPrintFlows || DoPrintBins)
          if(MPDEBUG) cout<<"myPlotInteractive 19 "<<endl;
          // END OF 1D HIST
        } else if(varDim==2) { 
          cout << "Plotting Histogram: " << ivar+jcat << "\t" << histnamebase << endl;
          outname = "hist_"; outname += plotvarnames[ivar];//histnamebase;
          TH2F * dummyhist = new TH2F("dummy","",1,-9999,-9998,1,-9999,-9998);
          dummyhist->Fill(-9998.5,-9998.5);
          dummyhist->SetStats(0);
          if(xaxislabels[ivar] != "xaxis")dummyhist->GetXaxis()->SetTitle(xaxislabels[ivar]);
          if(yaxislabels[ivar] != "yaxis")dummyhist->GetYaxis()->SetTitle(yaxislabels[ivar]);
          dummyhist->GetXaxis()->SetTitleSize(0.06);
          dummyhist->GetXaxis()->SetTitleOffset(0.75);
          dummyhist->GetYaxis()->SetTitleSize(0.06);
          dummyhist->GetYaxis()->SetTitleOffset(0.75);
          for(Int_t iInd=0;iInd!=NIND;++iInd){
            hist2D[iInd] = new TH2F();
            if(omitInd[iInd]%2==1)continue;
            *(hist2D[iInd]) = *((TH2F*)gROOT->FindObject(myHistname[iInd])); 
            hist2D[iInd]->GetXaxis()->SetLabelSize(0.04);
            hist2D[iInd]->GetYaxis()->SetLabelSize(0.04);
            if(xaxislabels[ivar] != "xaxis")hist2D[iInd]->GetXaxis()->SetTitle(xaxislabels[ivar]);
            if(yaxislabels[ivar] != "yaxis")hist2D[iInd]->GetYaxis()->SetTitle(yaxislabels[ivar]);
            hist2D[iInd]->GetXaxis()->SetTitleSize(0.06);
            hist2D[iInd]->GetXaxis()->SetTitleOffset(0.75);
            hist2D[iInd]->GetYaxis()->SetTitleSize(0.06);
            hist2D[iInd]->GetYaxis()->SetTitleOffset(0.75);
          }

          if(DoPrintFlows || DoPrintBins) {
            for(Int_t iInd=0;iInd!=NIND;++iInd){
              if(omitInd[iInd]%2==1)continue;
              Int_t NBinsX = hist2D[iInd]->GetNbinsX();
              Int_t NBinsY = hist2D[iInd]->GetNbinsY();
              float underflowX = hist2D[iInd]->Integral(0,0,0,NBinsY+1);
              float overflowX = hist2D[iInd]->Integral(NBinsX+1,NBinsX+1,0,NBinsY+1);
              float underflowY = hist2D[iInd]->Integral(0,NBinsX+1,0,0);
              float overflowY = hist2D[iInd]->Integral(0,NBinsX+1,NBinsY+1,NBinsY+1);
              float integral=hist2D[iInd]->Integral(0,NBinsX+1,0,NBinsY+1);
              if(underflowX>0 || overflowX>0 || underflowY>0 || overflowY>0) {
                std::cout << "\n" << myLabel[iInd] << std::endl;
                if(underflowX>0)std::cout << "\tunderflow X\t" << setw(9) << underflowX << "\t("<<underflowX/integral<<")" << std::endl;
                if(overflowX>0)std::cout << "\toverflow X\t" << setw(9) << overflowX << "\t("<<overflowX/integral<<")" << std::endl;
                if(underflowY>0)std::cout << "\tunderflow Y\t" << setw(9) << underflowY << "\t("<<underflowY/integral<<")" << std::endl;
                if(overflowY>0)std::cout << "\toverflow Y\t" << setw(9) << overflowY << "\t("<<overflowY/integral<<")" << std::endl;
              }
              //std::cout << myLabel[iInd] << "\ttotal\t" << integral << std::endl;
            }
          }//end if(DoPrintFlows || DoPrintBins)

          bool didfirst=false;
          for(Int_t iInd=0;iInd!=NIND;++iInd) {
            if(omitInd[myStackOrder[iInd]]%2==1)continue;
            hist2D[myStackOrder[iInd]]->SetTitle(histnamebase);
            hist2D[myStackOrder[iInd]]->SetStats(DoStats);
            hist2D[myStackOrder[iInd]]->SetMarkerColor(myColor[myStackOrder[iInd]]);
            hist2D[myStackOrder[iInd]]->SetMarkerStyle(myMarker[myStackOrder[iInd]]);
            //FIXME if(itypePI[infoind[iInd]]==0 && myMarker[myStackOrder[iInd]]>19 && myMarker[myStackOrder[iInd]]<24) hist2D[myStackOrder[iInd]]->SetMarkerStyle(myMarker[myStackOrder[iInd]]+4); 
            if(itypePI[iInd]==0 && myMarker[myStackOrder[iInd]]>19 && myMarker[myStackOrder[iInd]]<24) hist2D[myStackOrder[iInd]]->SetMarkerStyle(myMarker[myStackOrder[iInd]]+4); 
            if(didfirst) {
	      cout << "9911 drawing here 1" << endl; 
              hist2D[myStackOrder[iInd]]->Draw("same");
	    }
            else {
	      cout << "9912 drawing here 1" << endl; 
              hist2D[myStackOrder[iInd]]->Draw();
	    }
	    didfirst=true;
          }
          //dummyhist->Draw("same");
          //add hist to legend
          for(Int_t iInd=0;iInd!=NIND;++iInd){
            Int_t reverseInd = myStackOrder[NIND-1-iInd];// makes top of legend = top of stack
            if(omitInd[reverseInd]%2==1)continue;
            if(omitInd[reverseInd]==2)continue;
            legend->AddEntry(hist2D[reverseInd],myLabeldisplay[reverseInd],"p");
          }
          // END OF 2D HIST
        } else if(varDim==3) { 
          cout << "Plotting Profile Histogram: " << ivar+jcat << "\t" << histnamebase << endl;
          outname = "hist_"; outname += plotvarnames[ivar];//histnamebase;
          for(Int_t iInd=0;iInd!=NIND;++iInd) {
            histProf[iInd] = new TProfile();
            //FIXME if(omitInd[iInd]%2==1 || itypePI[infoind[iInd]]<0)continue;
            if(omitInd[iInd]%2==1 || itypePI[iInd]>0)continue;
            *(histProf[iInd]) = *((TProfile*)gROOT->FindObject(myHistname[iInd])); 
            histProf[iInd]->GetXaxis()->SetLabelSize(0.04);
            histProf[iInd]->GetYaxis()->SetLabelSize(0.04);
            if(xaxislabels[ivar] != "xaxis")histProf[iInd]->GetXaxis()->SetTitle(xaxislabels[ivar]);
            if(yaxislabels[ivar] != "yaxis")histProf[iInd]->GetYaxis()->SetTitle(yaxislabels[ivar]);
            histProf[iInd]->GetXaxis()->SetTitleSize(0.06);
            histProf[iInd]->GetXaxis()->SetTitleOffset(0.75);
            histProf[iInd]->GetYaxis()->SetTitleSize(0.06);
            histProf[iInd]->GetYaxis()->SetTitleOffset(0.75);
          }

          bool didfirst=false;
          for(Int_t iInd=NIND-1;iInd>=0;--iInd){
            //FIXME if(omitInd[myStackOrder[iInd]]%2==1 || itypePI[infoind[iInd]]<0)continue;
            if(omitInd[myStackOrder[iInd]]%2==1 || itypePI[iInd]>0)continue;
            histProf[myStackOrder[iInd]]->SetTitle(histnamebase);
            histProf[myStackOrder[iInd]]->SetStats(DoStats);
            histProf[myStackOrder[iInd]]->SetMarkerColor(myColor[myStackOrder[iInd]]);
            histProf[myStackOrder[iInd]]->SetMarkerStyle(27);
            //FIXME if(itypePI[infoind[iInd]]==0) histProf[myStackOrder[iInd]]->SetMarkerStyle(20);
            if(itypePI[iInd]==0) histProf[myStackOrder[iInd]]->SetMarkerStyle(20);
            //if(itypePI[infoind[iInd]]==0 && myMarker[myStackOrder[iInd]]>19 && myMarker[myStackOrder[iInd]]<24) histProf[myStackOrder[iInd]]->SetMarkerStyle(myMarker[myStackOrder[iInd]]+4); 
            histProf[myStackOrder[iInd]]->SetLineStyle(1);
            //FIXME if(itypePI[infoind[iInd]]==0) histProf[myStackOrder[iInd]]->SetLineWidth(2);
            if(itypePI[iInd]==0) histProf[myStackOrder[iInd]]->SetLineWidth(2);
            if(didfirst) {
	      cout << "9913 drawing here 1" << endl; 
              histProf[myStackOrder[iInd]]->Draw("same");
	    }
            else {
	      cout << "9914 drawing here 1" << endl; 
              histProf[myStackOrder[iInd]]->Draw();
	    }
            didfirst=true;
          }
          for(Int_t iInd=NIND-1;iInd>=0;--iInd){
            Int_t reverseInd = myStackOrder[NIND-1-iInd];// makes top of legend = top of stack
            //FIXME if(omitInd[reverseInd]%2==1 || itypePI[infoind[iInd]]<0)continue;
            if(omitInd[reverseInd]%2==1 || itypePI[iInd]>0)continue;
            if(omitInd[reverseInd]==2)continue;
            legend->AddEntry(histProf[reverseInd],myLabeldisplay[reverseInd],"p");
          }
          // END OF Profile HIST
        }

        //draw legend
          if(MPDEBUG) cout<<"myPlotInteractive 19.1 "<<endl;
        if(DoLegend) legend->Draw();
          if(MPDEBUG) cout<<"myPlotInteractive 19.2 "<<endl;
          if(MPDEBUG) cout<<"DoStack " << DoStack <<endl;
        if(DoStack) {
            for(int iline=0;iline!=NStackLines;++iline) {
              if(DoPopSig && histStackLineIsSig[iline])continue;
              histStackLine[iline].SetLineColor(kBlack);
              histStackLine[iline].SetLineWidth(lineWidth);
              histStackLine[iline].SetFillStyle(0);
              
              histStackLine[iline].Draw("hist same");
            }
            // plot signal lines last so background lines won't be in front
            for(int iline=0;iline!=NStackLines;++iline) {
              if( !(DoPopSig && histStackLineIsSig[iline]))continue;
	      //MARCOMARCO
              //histStackLine[iline].SetFillStyle(1001);
              histStackLine[iline].SetFillStyle(0);
	      //MARCOMARCO
              histStackLine[iline].SetLineColor(histStackLineColor[iline]);
              //histStackLine[iline].SetLineColor(kRed);
              histStackLine[iline].SetLineWidth(lineWidth);
              histStackLine[iline].Draw("hist same");
            }

	    //MARCOMARCO
            if(DoData) {
              if(myStack->GetHists() || mySigStack->GetHists()){
            cout << "9920 drawing here 1" << endl; 
                myDataStack->Draw("epsame");
              } 
	    }
        }
	float thishistmax =0.;
	//MARCO FIXME
	//float thishistmax = std::max(myStack->GetMaximum(),mySigStack->GetMaximum());



          if(MPDEBUG) cout<<"myPlotInteractive 19.3 "<<endl;
   if(sethistmax>0) thishistmax = sethistmax;
   TLine *cutline = new TLine();
   cutline->SetLineWidth(5);
   cutline->SetLineColor(632); //red
   cutline->DrawLine(linex,0.0,linex,thishistmax);
          if(MPDEBUG) cout<<"myPlotInteractive 19.4 "<<endl;

   if(textx2 <= textx1 || texty2 <= texty1) {
     textx1 = 0.325;
     textx2 = 0.6;
     texty1 = 0.83;
     texty2 = 0.9;
   }
          if(MPDEBUG) cout<<"myPlotInteractive 19.5 "<<endl;
   TPaveText *cmspavetext = new TPaveText(textx1,texty1,textx2,texty2,"brNDC");
   cmspavetext->SetTextFont(62);
   cmspavetext->SetTextSize(0.0425); //MARCOMARCO was 0.5
   cmspavetext->SetBorderSize(0);
   cmspavetext->SetLineColor(0);
   cmspavetext->SetLineStyle(0);
   cmspavetext->SetLineWidth(0);
   cmspavetext->SetFillColor(0);
   cmspavetext->SetFillStyle(0);
          if(MPDEBUG) cout<<"myPlotInteractive 19.6 "<<endl;
   TText *cmstext = cmspavetext->AddText("CMS Preliminary 2011");
   cmstext = cmspavetext->AddText("\n\n\n\n\n\n\n\n");
   //cmstext = cmspavetext->AddText("");
   //cmstext = cmspavetext->AddText("");
   cmstext = cmspavetext->AddText("#sqrt{s}=7 TeV L_{int} = 4.69 fb^{-1}");
   cmspavetext->Draw();
          if(MPDEBUG) cout<<"myPlotInteractive 19.7 "<<endl;
        //update canvas
        ch->Update();
          if(MPDEBUG) cout<<"myPlotInteractive 19.8 "<<endl;
        if(DoSnB) {
          ch2->Update();
          ch3->Update();
        }
        if(DoSvB) {
          ch4->Update();
          ch5->Update();
        }
        }// end for-loop (jcat)
    }// END if(DoRedraw)
    DoRedraw=true;
    ipair = myGetInput(ipair);
    ivar = ipair.first;
    icat = ipair.second;


  } 

  while(0 <= ivar && ivar < Nvar);
  //********************************************
  
}


std::pair<int,int> LoopAll::myGetInput(std::pair<int,int> ipair) {
//  for(int iiii=0;iiii!=Nvar;++iiii) {
//    std::cout << "doplot["<<iiii<<"]: " << doplot[iiii] << std::endl;
//  }
  int ivar=ipair.first;
  int icat=ipair.second;

  
  int ncat=histoncat[ivar];
  std::pair<int,int> newPrevPair(ivar,icat);
  TString sIn;
  Int_t iIn;
  TString logstring="";
  if(DoLog)logstring="_log";
  TString alpha[]={"aa","ab","ac","ad","ae","af","ag","ah","ai","aj","ak","al","am","an","ao","ap","aq","ar","as","at","au","av","aw","ax","ay","az",
		   "ba","bb","bc","bd","be","bf","bg","bh","bi","bj","bk","bl","bm","bn","bo","bp","bq","br","bs","bt","bu","bv","bw","bx","by","bz",
		   "ca","cb","cc","cd","ce","cf","cg","ch","ci","cj","ck","cl","cm","cn","co","cp","cq","cr","cs","ct","cu","cv","cw","cx","cy","cz",
		   "da","db","dc","dd","de","df","dg","dh","di","dj","dk","dl","dm","dn","do","dp","dq","dr","ds","dt","du","dv","dw","dx","dy","dz",
		   "ea","eb","ec","ed","ee","ef","eg","eh","ei","ej","ek","el","em","en","eo","ep","eq","er","es","et","eu","ev","ew","ex","ey","ez",
		   "fa","fb","fc","fd","fe","ff","fg","fh","fi","fj","fk","fl","fm","fn","fo","fp","fq","fr","fs","ft","fu","fv","fw","fx","fy","fz",
		   "ga","gb","gc","gd","ge","gf","gg","gh","gi","gj","gk","gl","gm","gn","go","gp","gq","gr","gs","gt","gu","gv","gw","gx","gy","gz"};
  if(DoWriteAll) {
    //ch->Print("allhists/"+outname+logstring+"."+allext);
    if(doplot[ivar]>0) {
      ch->Print("allhists/"+alpha[ivar]+outname+logstring+"."+allext);
    }
    while(doplot[++ivar]<=0 && ivar < Nvar)  {std::cout<<"not writing: "<<plotvarnames[ivar]<<std::endl;}
        
    DoBigLegend=false;
  } else {
    cout << "Enter command (h for help): ";
    std::cin >> sIn;
    iIn = atoi(sIn);
    if(iIn == 0) {
      // move through histograms
      if(sIn == "q"||sIn==".q")           {ivar = -1;}
      //else if(sIn == "n")						      {if(DoMulti){ivar++;}else{icat++;}}
      else if(sIn == "n")						      {if(DoMulti){while(doplot[++ivar]==0&&ivar<Nvar){std::cout<<"skipping: "<<plotvarnames[ivar]<<std::endl;}}else{icat++;}}
      else if(sIn=="j"||sIn=="jump")      {if(DoMulti){std::cin >> ivar;} else {std::cin >> ivar; std::cin >> icat;}}
      //else if(sIn == "p")						      {if(DoMulti){ivar--;}else{icat--;}}
      else if(sIn == "p")						      {if(DoMulti){while(doplot[--ivar]==0&&ivar<Nvar){std::cout<<"skipping: "<<plotvarnames[ivar]<<std::endl;}}else{icat--;}}
      else if(sIn == "b")				    		  {ivar=prevPair.first;icat=prevPair.second;}
        // legend
      else if(sIn == "legend")            {DoLegend=!DoLegend;}
      else if(sIn == "bigleg")            {DoBigLegend=!DoBigLegend;
        if(DoMulti && DoBigLegend){if(nrows*ncolumns>1){std::cout << "which hist? (1-" << nrows*ncolumns << "): ";std::cin >> BigLegendInd;}}
        if(BigLegendInd < 1 || BigLegendInd > nrows*ncolumns)BigLegendInd=1; }
      else if(sIn == "tr")                {DoLegendT=true;DoLegendR=true;}
      else if(sIn == "tl")                {DoLegendT=true;DoLegendR=false;}
      else if(sIn == "br")                {DoLegendT=false;DoLegendR=true;}
      else if(sIn == "bl")                {DoLegendT=false;DoLegendR=false;}
      // manipulate histograms
      else if(sIn == "stack")             {DoStack=!DoStack;DoFill=DoStack;}
      else if(sIn == "popsig")            {if(DoStack)DoPopSig=!DoPopSig;}
      else if(sIn == "data")              {DoData=!DoData;}
      else if(sIn == "scale")             {DoScale=!DoScale;}
      else if(sIn == "scalesig")          {DoScale=true;std::cin>>sigScale; }
      else if(sIn == "int")               {DoInt=!DoInt;DoRevInt=false;DoOverFlow=true;DoUnderFlow=true;}
      else if(sIn == "rint")              {DoRevInt=!DoRevInt;DoInt=false;DoOverFlow=true;DoUnderFlow=true;}
      else if(sIn == "svn")               {DoSvN=!DoSvN;DoInt = (!DoSvN||DoRevInt||DoInt)?DoInt:true; DoOverFlow=true;DoUnderFlow=true;}
      else if(sIn == "uflow")             {DoUnderFlow=!DoUnderFlow;}
      else if(sIn == "oflow")             {DoOverFlow=!DoOverFlow;}
      else if(sIn == "uof")					      {DoOverFlow=!DoOverFlow;DoUnderFlow=!DoUnderFlow;}
      else if(sIn == "rebin")             {DoRebin=!DoRebin;if(DoRebin)std::cin >> NReBin;}
      else if(sIn == "setmax")            {std::cout << "Enter Max: ";std::cin >> sethistmax;}
      else if(sIn == "setmin")            {std::cout << "Enter Min: ";std::cin >> sethistmin;}
      else if(sIn == "placeline")         {std::cout << "Enter x1 (old "<<linex<<"): ";std::cin >> linex;}
      else if(sIn == "placetext")         {std::cout << "Enter x1,x2,y1,y2 (old "<<textx1<<","<<textx2<<","<<texty1<<","<<texty2<<"): ";std::cin >> textx1 >> textx2 >> texty1 >> texty2;}
      else if(sIn == "placelegend" || sIn == "placeleg")  {std::cout << "Enter x1,x2,y1,y2 (old "<<legendx1<<","<<legendx2<<","<<legendy1<<","<<legendy2<<"): ";std::cin >> legendx1 >> legendx2 >> legendy1 >> legendy2;}
      else if(sIn == "stats")             {DoStats=!DoStats;}
      else if(sIn == "sumw2")             {DoSumw2=!DoSumw2;}
      else if(sIn == "log")               {DoLog=!DoLog;}
      else if(sIn == "fill")              {DoFill=!DoFill;}
      else if(sIn == "fillsig")           {DoFillSig=!DoFillSig;}
      else if(sIn == "rmtitles")          {RmTitles=!RmTitles;}
      else if(sIn == "cats")              {DoCats=!DoCats;DoMulti=DoCats; icat=0;}
      else if(sIn == "multi")             {if(!DoCats)DoMulti=!DoMulti; DoCats=false;icat=0; if(DoMulti){std::cin >> nrows;std::cin >> ncolumns;} else{nrows=1;ncolumns=1;}}
      else if(sIn == "startcat")          {cin >> StartCat;}
      else if(sIn == "sb")                {DoSnB=!DoSnB;
        if(DoSnB) {
          if(DoSvB) { 
            DoSvB = false;
            ch4->Close();
            ch5->Close();
          }
          ch2 = new TCanvas("ch2","ch2",10,ySize+55,xSize,ySize);
          ch3 = new TCanvas("ch3","ch3",xSize+20,10,xSize,ySize);
        } else {
          ch2->Close();
          ch3->Close();
        }
      }
      else if(sIn == "svb")                {DoSvB=!DoSvB;
        if(DoSvB) {
          if(DoSnB) { 
            DoSnB = false;
            ch2->Close();
            ch3->Close();
          }
          ch4 = new TCanvas("ch4","ch4",xSize+20,10,xSize,ySize); 
          ch5 = new TCanvas("ch5","ch5",10,ySize+55,xSize,ySize); 
        } else {
          ch4->Close();
          ch5->Close();
        }
      }
      // manipulate canvas
      else if(sIn == "grid")              {DoGridX=!DoGridX;DoGridY=!DoGridY;}
      else if(sIn == "gridx")             {DoGridX=!DoGridX;}
      else if(sIn == "gridy")             {DoGridY=!DoGridY;}
      else if(sIn == "fitsig")            {DoFitSig=!DoFitSig;if(DoFitSig){DoFit=true;std::cin >> fitXmin;std::cin >> fitXmax;}}
      else if(sIn == "fit")             {DoFit=!DoFit;}
      else if(sIn == "resize")            {std::cin >> xSize;std::cin >> ySize;}
      // other
      else if(sIn == "omit")              {std::cout<<"file\tomit\tSample"<<std::endl;
        for(Int_t i=0;i!=NIND;++i) { std::cout<<i<<"\t"<<omitInd[i]<<"\t"<<mySampleString[i]<<std::endl; }
        Int_t ifile; std::cout<<"which files? (-1 to end): ";
        do{ 
          std::cin >> ifile;
          if(ifile>-1 && ifile < NIND) {
            if(omitInd[ifile]%2==0) {
              omitInd[ifile]=1; 
            } else  {
              if(omitdefault[ifile]%2==0) {
                omitInd[ifile]=omitdefault[ifile]; 
              } else {
                omitInd[ifile]=0;
              }
            }
          }
        }while(ifile>-1); 
        } else if(sIn == "setsig")            {for(Int_t i=0;i!=NIND;++i)mySgnl[i]=false;std::cout<<"file\tSample"<<std::endl;
        for(Int_t i=0;i!=NIND;++i) {if(allSgnl[i])std::cout<<i<<"\t"<<mySampleString[i]<<std::endl; }
        Int_t ifile; std::cout<<"which files? (-1 to end): ";do{std::cin >> ifile;if(ifile>-1 && ifile < NIND)mySgnl[ifile]=true; }while(ifile>-1); }
      else if(sIn == "setback")            {for(Int_t i=0;i!=NIND;++i)myBkgd[i]=false;std::cout<<"file\tSample"<<std::endl;
        for(Int_t i=0;i!=NIND;++i) {if(allBkgd[i])std::cout<<i<<"\t"<<mySampleString[i]<<std::endl; }
        Int_t ifile; std::cout<<"which files? (-1 to end): ";do{std::cin >> ifile;if(ifile>-1 && ifile < NIND)myBkgd[ifile]=true; }while(ifile>-1); }
      else if(sIn == "write")             {DoRedraw=false;
        std::cout << "enter extension: ";TString ext;std::cin >> ext;ch->Print(outname+logstring+"."+ext);}
      //else if(sIn == "writeall")          {DoWriteAll=true;/////DoCats=true;////DoMulti=true;icat=0;ivar=0;DoBigLegend=true;BigLegendInd=1;
      else if(sIn == "writeall")          {DoWriteAll=true;DoMulti=true;icat=0;ivar=0;DoBigLegend=true;BigLegendInd=1;
        std::cout << "enter extension: ";std::cin >> allext;}
      else if(sIn == "writesob")          {DoRedraw=false;
        if(DoSnB) {std::cout << "enter extension: ";TString ext;std::cin >> ext;ch3->Print(outnamesob+logstring+"."+ext);}
        else{std::cout << "nowrite - not doing s/b" << std::endl;}}
      else if(sIn == "ls")					      {DoRedraw=false;for(int iv=0;iv!=Nvar;++iv)std::cout << "varname["<<iv<<"]: "<<plotvarnames[iv]<<std::endl;}
      else if(sIn == "pf" || sIn == "printflows")				{DoPrintFlows=!DoPrintFlows;}
      else if(sIn == "pb" || sIn == "printbins")					{DoPrintBins=!DoPrintBins;}
      else if(sIn == "h" || sIn == "help")    {DoRedraw=false;cout << OPTIONS << endl;}
    } else { ivar += iIn; }

    prevPair=newPrevPair;
    if(!DoMulti) {
      if(icat < 0) {
        icat=0;
        int newncat=ncat;
        if(ivar>0) newncat=histoncat[ivar-1];
        ivar--;
        if(newncat>0)icat=newncat-1;
      } else if(icat>=ncat && icat>0) {
        icat=0;
        ivar++;
      }
    }

  }
  

  std::pair<int,int> returnpair(ivar,icat);
  return returnpair;
}

void LoopAll::myPlotInteractiveSetup(TString hsmallname, TString tag) {
  cout<<" entring myPlotInteractiveSetup "<<hsmallname<<"  "<<tag<<endl;
  TFile *file1;
  cout<<"open file"<<endl;
    
  if(hsmallname!="") {
    if(MPDEBUG)std::cout << "PlotInteractive DEBUG 01 " << std::endl;
    file1 = TFile::Open(hsmallname);
    file1->cd();

    if(MPDEBUG)std::cout << "PlotInteractive DEBUG 02 " << std::endl;
    cout<<"get plotvariables"<<endl;
    plotvariables = (TTree*)file1->Get("plotvariables");
    cout<<"get inputfiles"<<endl;
    inputfiles = (TTree*)file1->Get("inputfiles");

    // plotvariables information BEGIN
    cout<<" set branch addresses"<<endl;
    plotvariables->SetBranchAddress("Nvar",&Nvar,&b_Nvar);
    plotvariables->SetBranchAddress("typplotall",&typplotall,&b_typplotall);
    plotvariables->SetBranchAddress("doplot",&doplot,&b_doplot);
    plotvariables->SetBranchAddress("h2d",&h2d,&b_h2d);
    plotvariables->SetBranchAddress("typplot",&typplot,&b_typplot);
    plotvariables->SetBranchAddress("histoncat",&histoncat,&b_histoncat);
    plotvariables->SetBranchAddress("histoncatindtonames",&histoncatindtonames,&b_histoncatindtonames);
    plotvariables->SetBranchAddress("nbinsx",&nbinsx,&b_nbinsx);
    plotvariables->SetBranchAddress("nbinsy",&nbinsy,&b_nbinsy);
    plotvariables->SetBranchAddress("lowlim",&lowlim,&b_lowlim);
    plotvariables->SetBranchAddress("highlim",&highlim,&b_highlim);
    plotvariables->SetBranchAddress("lowlim2",&lowlim2,&b_lowlim2);
    plotvariables->SetBranchAddress("highlim2",&highlim2,&b_highlim2);
    plotvariables->SetBranchAddress("xaxislabels",&tca_xaxislabels,&b_xaxislabels);
    plotvariables->SetBranchAddress("yaxislabels",&tca_yaxislabels,&b_yaxislabels);
    plotvariables->SetBranchAddress("plotvarnames",&tca_plotvarnames,&b_plotvarnames);
    plotvariables->SetBranchAddress("Nvarcats",&Nvarcats,&b_Nvarcats);
    plotvariables->SetBranchAddress("catid",&catid,&b_catid);
    plotvariables->SetBranchAddress("ncats",&ncats,&b_ncats);
    plotvariables->SetBranchAddress("plotvarcatnames",&tca_plotvarcatnames,&b_plotvarcatnames);

    cout<<"get branches"<<endl;
    b_Nvar = plotvariables->GetBranch("Nvar");
    b_typplotall = plotvariables->GetBranch("typplotall");
    b_doplot = plotvariables->GetBranch("doplot");
    b_h2d = plotvariables->GetBranch("h2d");
    b_typplot = plotvariables->GetBranch("typplot");
    b_histoncat = plotvariables->GetBranch("histoncat");
    b_histoncatindtonames = plotvariables->GetBranch("histoncatindtonames");
    b_nbinsx = plotvariables->GetBranch("nbinsx");
    b_nbinsy = plotvariables->GetBranch("nbinsy");
    b_lowlim = plotvariables->GetBranch("lowlim");
    b_highlim = plotvariables->GetBranch("highlim");
    b_lowlim2 = plotvariables->GetBranch("lowlim2");
    b_highlim2 = plotvariables->GetBranch("highlim2");
    b_xaxislabels = plotvariables->GetBranch("xaxislabels");
    b_yaxislabels = plotvariables->GetBranch("yaxislabels");
    b_plotvarnames = plotvariables->GetBranch("plotvarnames");
    //b_Nvarcats = plotvariables->GetBranch("Nvarcats");
    //b_catid = plotvariables->GetBranch("catid");
    //b_ncats = plotvariables->GetBranch("ncats");
    //b_plotvarcatnames = plotvariables->GetBranch("plotvarcatnames");

    cout<<"get Entries"<<endl;
    b_Nvar->GetEntry();
    b_typplotall->GetEntry();
    b_doplot->GetEntry();
    b_h2d->GetEntry();
    b_typplot->GetEntry();
    b_histoncat->GetEntry();
    b_histoncatindtonames->GetEntry();

    b_nbinsx->GetEntry();
    b_nbinsy->GetEntry();
    b_lowlim->GetEntry();
    b_highlim->GetEntry();
    b_lowlim2->GetEntry();
    b_highlim2->GetEntry();
    b_xaxislabels->GetEntry();
    b_yaxislabels->GetEntry();
    b_plotvarnames->GetEntry();
    //b_Nvarcats->GetEntry();
    //b_catid->GetEntry();
    //b_ncats->GetEntry();
    //b_plotvarcatnames->GetEntry();

    // plotvariables information END
    if(MPDEBUG)std::cout << "PlotInteractive DEBUG 03 " << std::endl;
    // inputfiles information BEGIN
    cout<<"set branches"<<endl;
    inputfiles->SetBranchAddress("nfiles",&nfiles,&b_nfiles);
    inputfiles->SetBranchAddress("nindfiles",&nindfiles,&b_nindfiles);
    inputfiles->SetBranchAddress("intlumi",&intlumi,&b_intlumi);
    inputfiles->SetBranchAddress("makeOutputTree",&makeOutputTree,&b_makeOutputTree);
    inputfiles->SetBranchAddress("histfilename",&tca_histfilename,&b_histfilename);
    inputfiles->SetBranchAddress("itype",&itypePI,&b_itype);
    inputfiles->SetBranchAddress("histoind",&histoind,&b_histoind);
    inputfiles->SetBranchAddress("infoind",&infoind,&b_infoind);
    inputfiles->SetBranchAddress("histoplotit",&histoplotit,&b_histoplotit);
    inputfiles->SetBranchAddress("ntot",&ntot,&b_ntot);
    inputfiles->SetBranchAddress("nred",&nred,&b_nred);
    inputfiles->SetBranchAddress("lumi",&lumi,&b_lumi);
    inputfiles->SetBranchAddress("xsec",&xsec,&b_xsec);
    inputfiles->SetBranchAddress("kfactor",&kfactor,&b_kfactor);
    inputfiles->SetBranchAddress("scale",&scale,&b_scale);
    inputfiles->SetBranchAddress("inshortnames",&tca_inshortnames,&b_inshortnames);
    inputfiles->SetBranchAddress("infilenames",&tca_infilenames,&b_infilenames);

    cout<<"get branches"<<endl;
    b_nfiles = inputfiles->GetBranch("nfiles");
    b_nindfiles = inputfiles->GetBranch("nindfiles");
    b_intlumi = inputfiles->GetBranch("intlumi");
    b_makeOutputTree = inputfiles->GetBranch("makeOutputTree");
    b_histfilename = inputfiles->GetBranch("histfilename");
    b_itype = inputfiles->GetBranch("itype");
    b_histoind = inputfiles->GetBranch("histoind");
    b_infoind = inputfiles->GetBranch("infoind");
    b_histoplotit = inputfiles->GetBranch("histoplotit");
    b_ntot = inputfiles->GetBranch("ntot");
    b_nred = inputfiles->GetBranch("nred");
    b_lumi = inputfiles->GetBranch("lumi");
    b_xsec = inputfiles->GetBranch("xsec");
    b_kfactor = inputfiles->GetBranch("kfactor");
    b_scale = inputfiles->GetBranch("scale");
    b_inshortnames = inputfiles->GetBranch("inshortnames");
    b_infilenames = inputfiles->GetBranch("infilenames");

    cout<<"get Entries"<<endl;
    b_nfiles->GetEntry();
    b_nindfiles->GetEntry();
    b_intlumi->GetEntry();
    b_makeOutputTree->GetEntry();
    b_histfilename->GetEntry();
    b_itype->GetEntry();
    b_histoind->GetEntry();
    b_infoind->GetEntry();
    b_histoplotit->GetEntry();
    b_ntot->GetEntry();
    b_nred->GetEntry();
    b_lumi->GetEntry();
    b_xsec->GetEntry();
    b_kfactor->GetEntry();
    b_scale->GetEntry();
    b_inshortnames->GetEntry();
    b_infilenames->GetEntry();

    if(MPDEBUG)std::cout << "nfiles: " << nfiles <<" " <<nindfiles<<std::endl;
    std::cout << "nfiles: " << nfiles <<" " <<nindfiles<<" "<<ntot[0]<<" "<<std::endl;

    nfiles=nindfiles;

    histfilename = ((TObjString*)tca_histfilename->At(0))->GetString();


    infilenames[0]=TString("Data");
    cout<<infilenames[0]<<endl;
    cout<<"before"<<endl;
    inshortnames[0]=TString(infilenames[0]);
    cout<<inshortnames[0]<<endl;
    cout<<"after"<<endl;


      for(Int_t ifile=0;ifile!=nfiles;ifile++) {
	//MARCO??? FIXME???
	//infilenames[ifile] = ((TObjmpUtil * mp, String*)tca_infilenames->At(ifile))->GetString();
	infilenames[ifile] = ((TObjString*)tca_infilenames->At(ifile))->GetString();
	inshortnames[ifile] = ((TObjString*)tca_inshortnames->At(ifile))->GetString();
      }
    
      /*
	nfiles=10;
    
	infilenames[0]=TString("Data");
	cout<<infilenames[0]<<endl;
	infilenames[1]=TString("Box25");
	infilenames[2]=TString("Box250");
	infilenames[3]=TString("DiPhotonJets");
	infilenames[4]=TString("GJet20_pf");
	infilenames[5]=TString("QCD40_fp");
	infilenames[6]=TString("QCD40_ff");
	infilenames[7]=TString("DYJets");
	infilenames[8]=TString("vbf_H_gg_150_pu2011");
	infilenames[9]=TString("wz_H_gg_150_pu2011");
	//infilenames[8]=TString("VBF150");
	//infilenames[9]=TString("WZH150");
	cout<<infilenames[9]<<endl;
	inshortnames[0]=TString(infilenames[0]);
	cout<<inshortnames[0]<<endl;
	inshortnames[1]=TString(infilenames[1]);
	inshortnames[2]=TString(infilenames[2]);
	inshortnames[3]=TString(infilenames[3]);
	inshortnames[4]=TString(infilenames[4]);
	inshortnames[5]=TString(infilenames[5]);
	inshortnames[6]=TString(infilenames[6]);
	inshortnames[7]=TString(infilenames[7]);
	inshortnames[8]=TString(infilenames[8]);
	inshortnames[9]=TString(infilenames[9]);
	cout<<inshortnames[9]<<endl;
      */

    cout<<"done fillingfile"<<endl;

    // fill strings from plotvariables
    for( Int_t ivar=0;ivar!=Nvar;++ivar) {
      plotvarnames[ivar] = ((TObjString*)tca_plotvarnames->At(ivar))->GetString();
      std::cout << "plotvarnames["<<ivar<<"]: " << plotvarnames[ivar] << std::endl;
      xaxislabels[ivar] = ((TObjString*)tca_xaxislabels->At(ivar))->GetString();
      yaxislabels[ivar] = ((TObjString*)tca_yaxislabels->At(ivar))->GetString();
      xaxislabels[ivar].ReplaceAll("@"," ");
      yaxislabels[ivar].ReplaceAll("@"," ");
    }
    int ivarcats=0;
    for( Int_t ivarcat=0;ivarcat!=Nvarcats;++ivarcat) {
      for(int iicat=0;iicat!=ncats[ivarcat];iicat++) {
        plotvarcatnames[ivarcat][iicat] = ((TObjString*)tca_plotvarcatnames->At(ivarcats++))->GetString();
      }
    }
    cout<<"done filling plotvariables"<<endl;

  }
    if(MPDEBUG)std::cout << "PlotInteractive DEBUG 04 " << std::endl;
  // inputfiles information END
  //NFILES = mp->nfiles;
  //NIND = mp->nindfiles;
  NFILES=nfiles;
  NIND=nindfiles;

  NIND=nfiles;


  xSizeMAX=1500;      
  ySizeMAX=950;     
  xSize=1500;
  ySize=600;
  nrows=1;
  ncolumns=1;
  sigScale=0;

  DoCats=true;
  DoLegend=true;
  DoLegendT=true;
  DoLegendR=true;
  DoMulti=true;
  DoScale=true;
  DoRedraw=true;
  DoLog=false;
  DoRebin=false;
  DoOverFlow=false;
  DoUnderFlow=false;
  DoStats=false;
  DoSumw2=false;
  DoStack=false;
  DoPopSig=false;
  DoData=true;
  DoFill=false;
  DoFillSig=false;
  DoBigLegend=false;
  DoGridX=false;
  DoGridY=false;
  DoInt=false;
  DoRevInt=false;
  DoSvN=false;
  DoFitSig=false;
  DoFit=false;
  DoPrintFlows=false;
  DoPrintBins=false;
  DoSnB=false;
  DoSvB=false;


  //if(MPDEBUG) {
    //std::cout << "mp->NFILES: " << mp->nfiles << std::endl;
    //std::cout << "mp->NIND: " << mp->nindfiles << std::endl;
    if(MPDEBUG)std::cout << "PlotInteractive DEBUG 05 " << std::endl;
    if(MPDEBUG)std::cout << "NFILES: " << nfiles << std::endl;
    if(MPDEBUG)std::cout << "NIND: " << nindfiles << std::endl;
    for(Int_t iFile=0;iFile!=NFILES;++iFile) {
      /*
      if(MPDEBUG) std::cout << "iFile: " << iFile << 
        "\titype: " << itypePI[iFile] << 
        //"\tmp itype: " << mp->itypePI[iFile] << 
        "\thistoind: " << histoind[iFile] << 
        "\tmp histoind: " << mp->histoind[iFile] << 
        "\tmp histoindfromfiles: " << mp->histoindfromfiles[iFile] << 
        "\thistoind[infoind]: " << histoind[infoind[iFile]] << 
        "\tinfoind: " << infoind[iFile] << 
        //"\tmp infoind: " << mp->infoind[iFile] << 
        "\tshortname: " << inshortnames[iFile] << 
        //"\tmp shortname: " << mp->filesshortnam[iFile] <<
        std::endl;
      */
    }
    //}
    if(MPDEBUG)std::cout << "PlotInteractive DEBUG 06 " << std::endl;

  bool defaultsetup=true;
  if(tag != "") {// read anaysis configuration from file
    defaultsetup=false;
    FILE * PIfile;
    if(tag == "Hgg") { PIfile=fopen("Marco/plotinteractiveHgg.dat","r");
    } else if(tag == "sgs1") { PIfile=fopen("sgs1_scripts/plotinteractivesgs1.dat","r");
    } else if(tag == "Hww") { PIfile=fopen("Hww_scripts/plotinteractiveHww.dat","r"); 
    } else if(tag == "Hzz") { PIfile=fopen("Hzz_scripts/plotinteractiveHzz.dat","r");
    } else if(tag == "Mwl") { PIfile=fopen("../Mwl_scripts/plotinteractiveMwl.dat","r");
    } else if(tag == "Elizabeth") { PIfile=fopen("/home/users/edusinberre/ActiveProjects/globe415/Elizabeth/analyze/plotinteractiveElizabeth.dat","r");
    } else if(tag == "estar") { PIfile=fopen("plotinteractiveestar.dat","r");
    } else if(tag == "threeSCs") { PIfile=fopen("interactiveplot3SC.dat","r");
    } else if(tag == "WAna") { PIfile=fopen("interactiveplotWAna.dat","r");
    } else if(tag == "Zee") { PIfile=fopen("plotinteractiveZee.dat","r");
//SED_END_PLOTINTERACTIVE
    }

    if(MPDEBUG)std::cout << "PlotInteractive DEBUG 07 " << std::endl;
    Int_t xSizeREAD,ySizeREAD;
    Int_t DoStackREAD,DoPopSigREAD,DoFillSigREAD,DoOverFlowREAD,DoUnderFlowREAD;
    Int_t DoLogREAD,DoGridREAD,DoStatsREAD,DoSumW2READ,DoCatsREAD;
    Int_t DoLegendREAD,DoLegendTREAD,DoLegendRREAD,NColsREAD;
    char name[100];
    char displayname[100];
    Int_t ind[NIND],stackorder[NIND],omit[NIND],color[NIND],marker[NIND];
    Float_t scale[NIND];

    int dummy = fscanf(PIfile,"xSize=%d ySize=%d\n",&xSizeREAD, &ySizeREAD);
    dummy = fscanf(PIfile,"Stack=%d PopSig=%d FillSig=%d Over=%d Under=%d\n",&DoStackREAD,&DoPopSigREAD,&DoFillSigREAD,&DoOverFlowREAD,&DoUnderFlowREAD);
    dummy =fscanf(PIfile,"Log=%d Grid=%d Stats=%d SumW2=%d DoCats=%d\n",&DoLogREAD,&DoGridREAD,&DoStatsREAD,&DoSumW2READ,&DoCatsREAD);
    dummy =fscanf(PIfile,"Legend=%d Top=%d Right=%d NCols=%d\n",&DoLegendREAD,&DoLegendTREAD,&DoLegendRREAD,&NColsREAD);
    for( int iInd=0;iInd!=NIND;++iInd) {
      //dummy = fscanf(PIfile,"name=%s ind=%d order=%d omit=%d scale=%f color=%d marker=%d\n", name, &ind[iInd], &stackorder[iInd], &omitdefault[iInd], &scale[iInd], &color[iInd], &marker[iInd]);
      dummy = fscanf(PIfile,"name=%s display=%s ind=%d order=%d omit=%d scale=%f color=%d marker=%d\n", name, displayname, &ind[iInd], &stackorder[iInd], &omitdefault[iInd], &scale[iInd], &color[iInd], &marker[iInd]);
      myLabel[iInd]=name;
      myLabeldisplay[iInd]=displayname;



      if(ind[iInd]!=histoind[infoind[iInd]]) {

        std::cout << "\nPROBLEM! int plotinteractive***.dat file - index does not match from inputfiles: ind:" << ind[iInd] << " should be " << histoind[infoind[iInd]]<<" [infoind[iInd] "<<infoind[iInd];
        std::cout << "\t. . . NOT Using Defualt PlotInteractive Setup . . . " << std::endl;
			

		      //MARCO FIXME												
	histoind[infoind[iInd]]=ind[iInd];

	infoind[iInd]=ind[iInd];



																					      //MARCO FIXME
																							      //defaultsetup=true;																			      //break;
      }
    }
    NColsLegend = NColsREAD;
    if(MPDEBUG)std::cout << "PlotInteractive DEBUG 08 " << std::endl;
    if(!defaultsetup) {
      xSize=xSizeREAD;      
      ySize=ySizeREAD;      

      DoStack=(bool)DoStackREAD;
      DoPopSig=(bool)DoPopSigREAD;
      DoFill=DoStack;
      DoFillSig=(bool)DoFillSigREAD;
      DoOverFlow=(bool)DoOverFlowREAD;
      DoUnderFlow=(bool)DoUnderFlowREAD;

      DoLog=(bool)DoLogREAD;
      DoGridX=(bool)DoGridREAD;
      DoGridY=(bool)DoGridREAD;
      DoStats=(bool)DoStatsREAD;
      DoSumw2=(bool)DoSumW2READ;

      DoLegend=(bool)DoLegendREAD;
      DoLegendT=(bool)DoLegendTREAD;
      DoLegendR=(bool)DoLegendRREAD;


    if(MPDEBUG)std::cout << "PlotInteractive DEBUG 09 " << std::endl;
      for( int iInd=0;iInd!=NIND;++iInd) {
    if(MPDEBUG)std::cout << "PlotInteractive DEBUG 09 - 01" << std::endl;
        myStackOrder[iInd] = stackorder[iInd];
    if(MPDEBUG)std::cout << "PlotInteractive DEBUG 09 - 02" << std::endl;
    std::cout << iInd << std::endl;
    std::cout << infoind[iInd] << std::endl;
    std::cout << inshortnames[iInd] << std::endl;
        mySampleString[iInd]=inshortnames[iInd];
    if(MPDEBUG)std::cout << "PlotInteractive DEBUG 09 - 03" << std::endl;
        omitInd[iInd]=omitdefault[iInd];
    if(MPDEBUG)std::cout << "PlotInteractive DEBUG 09 - 04" << std::endl;
        myScale[iInd]=scale[iInd];
    if(MPDEBUG)std::cout << "PlotInteractive DEBUG 09 - 05" << std::endl;
        //myLabel[iInd]=inshortnames[infoind[iInd]];
    if(MPDEBUG)std::cout << "PlotInteractive DEBUG 09 - 06" << std::endl;
        myMarker[iInd]=marker[iInd];
    if(MPDEBUG)std::cout << "PlotInteractive DEBUG 09 - 07" << std::endl;
        myColor[iInd]=color[iInd];
    if(MPDEBUG)std::cout << "PlotInteractive DEBUG 09 - 08" << std::endl;

      }
    if(MPDEBUG)std::cout << "PlotInteractive DEBUG 10 " << std::endl;
    }
    fclose(PIfile);
  }// END read anaysis configuration from file
    if(MPDEBUG)std::cout << "PlotInteractive DEBUG 11 " << std::endl;
  if(tag=="" || defaultsetup) {

    for(Int_t iInd=0;iInd!=NIND;++iInd){ 

      myStackOrder[iInd]=iInd;
      mySampleString[iInd]=inshortnames[infoind[iInd]];
      omitInd[iInd]=0;
      myScale[iInd]=1.;
      myLabel[iInd]=inshortnames[infoind[iInd]];
      myColor[iInd]=kGreen;

      if(itypePI[infoind[iInd]]<0) {
        myColor[iInd]=kRed;
      } else if(itypePI[infoind[iInd]]>0) {
        myColor[iInd]=kBlue+iInd%4;
      }
      myMarker[iInd]=iInd;
    }
  }

    if(MPDEBUG)std::cout << "PlotInteractive DEBUG 12 " << std::endl;
  for(Int_t iInd=0;iInd!=NIND;++iInd){ 
    myLabeldisplay[iInd].ReplaceAll("@"," ");
    allSgnl[iInd]=false;
    mySgnl[iInd]=false;
    myBkgd[iInd]=false;
    allBkgd[iInd]=false;

    cout<<iInd<<" "<<infoind[iInd]<<endl;

    if(itypePI[infoind[iInd]]<0) {
      allSgnl[iInd]=true;
      mySgnl[iInd]=true;
    } else if(itypePI[infoind[iInd]]>0) {
      allBkgd[iInd]=true;
      myBkgd[iInd]=true;
    }
  }

    
    if(MPDEBUG)std::cout << "PlotInteractive DEBUG 13 " << std::endl;
}

TString LoopAll::num_to_string(int myint) {
  std::stringstream stream;
  stream << myint;
  return (TString)stream.str();
}

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//============================  End Interactive Plotting  =============================
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
