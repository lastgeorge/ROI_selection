double cal_rms1(TH1F *htemp){
  //calculate rms, this is to be used for threshold purpose
  float rms = 0, rms1 = 0,rms2 = 0;
    
  // new method to calculate RMS
  int min1 =0,max1=0;
  for (int i=0;i!=htemp->GetNbinsX();i++){
    
    if (htemp->GetBinContent(i+1)>max1)
      max1 = int(htemp->GetBinContent(i+1));
    if (htemp->GetBinContent(i+1)<min1)
      min1 = int(htemp->GetBinContent(i+1));
    
  }
  TH1F *h6 = new TH1F("h6","h6",int(max1-min1+1),min1,max1+1);
  for (int i=0;i!=htemp->GetNbinsX();i++){
    h6->Fill(int(htemp->GetBinContent(i+1)));
  }
  if (h6->GetSum()>0){
    //calculate 0.16, 0.84 percentile ...  
    double xq;
    xq = 0.16;
    double par[3];
    h6->GetQuantiles(1,&par[0],&xq);
    xq = 0.84;
    h6->GetQuantiles(1,&par[1],&xq);
    xq = 0.5;
    h6->GetQuantiles(1,&par[2],&xq);
    rms = sqrt((pow(par[1]-par[2],2)+pow(par[0]-par[2],2))/2.);
  }
    
    delete h6;
    
    return rms;
}

  void cal_rms(){
    // TFile *file = new TFile("2D_display_3493_821_0.root");
    // TH2F *hu = (TH2F*)file->Get("hu_raw");
    // TH2F *hv = (TH2F*)file->Get("hv_raw");
    // TH2F *hw = (TH2F*)file->Get("hw_raw");

    // TFile *file = new TFile("1D_all.root");
    // TH2F *hu = (TH2F*)file->Get("hdecon_u");
    // TH2F *hv = (TH2F*)file->Get("hdecon_v");
    // TH2F *hw = (TH2F*)file->Get("hdecon_w");

    TFile *file = new TFile("2D_lf_all.root");
    TH2F *hu = (TH2F*)file->Get("hdecon_u");
    TH2F *hv = (TH2F*)file->Get("hdecon_v");
    TH2F *hw = (TH2F*)file->Get("hdecon_w");

    Int_t ntime = hu->GetNbinsY();
    Int_t nwire;

    ofstream outfile("result.dat");

    Int_t rebin = 6;

    nwire = hu->GetNbinsX();
    for (Int_t i=0;i!=nwire;i++){
      TH1F *htemp = new TH1F("htemp","htemp",ntime,0,ntime);
      //      htemp->Reset();
      for (Int_t j=0;j!=ntime;j++){
	htemp->SetBinContent(j+1,hu->GetBinContent(i+1,j+1));
      }
      htemp->Rebin(rebin);
      outfile << i << " " << cal_rms1(htemp) << endl;
      delete htemp;
    }
    
    nwire = hv->GetNbinsX();
    for (Int_t i=0;i!=nwire;i++){
      TH1F *htemp = new TH1F("htemp","htemp",ntime,0,ntime);
      
      for (Int_t j=0;j!=ntime;j++){
	htemp->SetBinContent(j+1,hv->GetBinContent(i+1,j+1));
      }
      htemp->Rebin(rebin);
      outfile << i << " " << cal_rms1(htemp) << endl;
      delete htemp;
    }

    nwire = hw->GetNbinsX();
    for (Int_t i=0;i!=nwire;i++){
      TH1F *htemp = new TH1F("htemp","htemp",ntime,0,ntime);
      
      for (Int_t j=0;j!=ntime;j++){
	htemp->SetBinContent(j+1,hw->GetBinContent(i+1,j+1));
      }
      htemp->Rebin(rebin);
      outfile << i << " " << cal_rms1(htemp) << endl;
      delete htemp;
    }

  }

