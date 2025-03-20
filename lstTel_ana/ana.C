#define ana_cxx
#include "ana.h"
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TMath.h>

void ana::Loop()
{

  TH1D *h1_nPhot = new TH1D("h1_nPhot","nPhot",100,0.0,100000); 
  TH1D *h1_PosX  = new TH1D("h1_PosX" ,"PosX" ,1000,-25000,25000); 
  TH1D *h1_PosY  = new TH1D("h1_PosY" ,"PosY" ,1000,-25000,25000); 
  //TH1D *h1_PosZ  = new TH1D("h1_PosZ" ,"PosZ" ,1000,0,1000); 
  //TH1D *h1_time  = new TH1D("h1_time" ,"time" ,100000,10.0,12.0); 

  TH1D *h1_R  = new TH1D("h1_R" ,"R" ,1000,0,200); 

  
  TH2D *h2_PosX_vs_PosY  = new TH2D("h2_PosX_vs_PosY" ,"PosX vs PosY" ,400,-25000,25000,400,-25000,25000); 

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"nentries = "<<nentries<<endl;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    h1_nPhot->Fill(nPhot);
    for(Int_t i = 0;i<nPhot;i++){
      h1_PosX->Fill(PosX[i]);
      //if(jentry == 10)
	h2_PosX_vs_PosY->Fill(PosX[i], PosY[i]);
      //h1_R->Fill(TMath::Sqrt(PosX[i]*PosX[i] + PosY[i]*PosY[i]));
    }
  }

  TCanvas *c1 = new TCanvas("c1","canva",10,10,1500,1000);
  //h1_nPhot->Draw();
        h2_PosX_vs_PosY->Draw("ZCOLOR");
  //TH1D *h1_PosX  = new TH1D("h1_PosX" ,"PosX" ,1000,-200000,200000); 
  //TH1D *h1_PosY  = new TH1D("h1_PosY" ,"PosY" ,1000,-200000,200000); 

  //h2_PosX_vs_PosY->Draw("ZCOLOR");

  
  TCanvas *c2 = new TCanvas("c2","canva",10,10,1500,1000);
  h1_PosX->Draw();
  
}
