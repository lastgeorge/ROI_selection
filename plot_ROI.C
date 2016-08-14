void restore_baseline(TH1F *htemp){
  //correct baseline 
  double max = htemp->GetMaximum();
  double min = htemp->GetMinimum();
  int nbin_b = max - min;
  Int_t nbin = htemp->GetNbinsX();
  if (nbin_b ==0) nbin_b = 1;
  TH1F *h1 = new TH1F("h1","h1",nbin_b,min,max);
  for (int j=0;j!=nbin;j++){
    h1->Fill(htemp->GetBinContent(j+1));
  }
  float ped = h1->GetMaximumBin()*(max-min)/(nbin_b*1.) + min;
  float ave=0,ncount = 0;
  
  for (int j=0;j!=nbin;j++){
    if (fabs(htemp->GetBinContent(j+1)-ped)<400){
      ave +=htemp->GetBinContent(j+1);
      ncount ++;
      }
    }
    if (ncount==0) ncount=1;
    ave = ave/ncount;
    
    for (int j=0;j!=nbin;j++){
      double content = htemp->GetBinContent(j+1);
      content -= ave;
      htemp->SetBinContent(j+1,content);
    }
    delete h1;
}



void plot_ROI(){
  
  Int_t plane = 0;
  //  Int_t run = 3493;
  Int_t run = 5366;
  Int_t wire_bin_min; 
  Int_t wire_bin_max; 
  Int_t time_bin_min; 
  Int_t time_bin_max; 
  wire_bin_min = 1170; 
  wire_bin_max = 1300;

  // wire_bin_min = 3850-2400; 
  // wire_bin_max = 3930-2400;
  time_bin_min = 6200;
  time_bin_max = 7800;

  // wire_bin_min = 1220; 
  // wire_bin_max = 1350;
  // time_bin_min = 6300;
  // time_bin_max = 7500;
  
  TFile *file = new TFile(Form("%d_2D_lf_all_25.root",run));

  TH2F *hdecon;
  TH2F *hdecon1;
  TH1F *htemp;

  
  
  if (plane == 0){
    hdecon = (TH2F*)file->Get("hdecon_u");
  }else if (plane == 1){
    hdecon = (TH2F*)file->Get("hdecon_v");
  }else if (plane == 2){
    hdecon = (TH2F*)file->Get("hdecon_w");
  }

  

  Int_t nwire = wire_bin_max - wire_bin_min + 1;
  Int_t ntime = time_bin_max - time_bin_min + 1;
 
  hdecon1 = new TH2F("hdecon1","hdecon1",nwire,0,nwire,ntime,0,ntime);
  htemp = new TH1F("htemp","htemp",ntime,0,ntime);
  
  for (Int_t i=wire_bin_min; i!= wire_bin_max+1; i++){
    for (Int_t j=time_bin_min; j!=time_bin_max+1; j++){
      htemp->SetBinContent(j-time_bin_min+1,hdecon->GetBinContent(i+1,j+1));
    }
    restore_baseline(htemp);
    
    for (Int_t j=time_bin_min; j!=time_bin_max+1; j++){
      hdecon1->SetBinContent(i-wire_bin_min+1,j-time_bin_min+1,htemp->GetBinContent(j-time_bin_min+1));
    }
  }
  
  hdecon1->Draw("COLZ");
  
  Double_t factor1 = 6;

  TFile *file1 = new TFile(Form("ROI_2D_25_%d.root",run));
  TTree *T = (TTree*)file1->Get("T");
  Int_t s_begin,s_end,s_channel,s_plane;
  T->SetBranchAddress("begin",&s_begin);
  T->SetBranchAddress("end",&s_end);
  T->SetBranchAddress("channel",&s_channel);
  T->SetBranchAddress("plane",&s_plane);
  for (Int_t i=0;i!=T->GetEntries();i++){
    T->GetEntry(i);
    if (plane == s_plane && s_channel >= wire_bin_min && s_channel <= wire_bin_max){
      TLine *l1 = new TLine(s_channel - wire_bin_min+0.5,s_begin*factor1-time_bin_min+0.5,s_channel - wire_bin_min+0.5,s_end*factor1-time_bin_min-0.5);
      l1->Draw("same");
    }
  }

  // TFile *file2 = new TFile(Form("ROI_1D_%d.root",run));
  // TTree *T1 = (TTree*)file2->Get("T");
  // T1->SetBranchAddress("begin",&s_begin);
  // T1->SetBranchAddress("end",&s_end);
  // T1->SetBranchAddress("channel",&s_channel);
  // T1->SetBranchAddress("plane",&s_plane);
  // for (Int_t i=0;i!=T1->GetEntries();i++){
  //   T1->GetEntry(i);
  //   if (plane == s_plane && s_channel >= wire_bin_min && s_channel <= wire_bin_max){
  //     TLine *l1 = new TLine(s_channel - wire_bin_min+0.5,s_begin*factor1-time_bin_min+0.5,s_channel - wire_bin_min+0.5,s_end*factor1-time_bin_min-0.5);
  //     l1->Draw("same");
  //     l1->SetLineColor(6);
  //   }
  // }

}
