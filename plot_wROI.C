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

Double_t local_ave(TH1F *h1, Int_t bin, Int_t width){
  Double_t sum1 = 0;
  Double_t sum2 = 0;
  
  for (Int_t i=-width;i<width+1;i++){
    Int_t current_bin = bin + i;

    while (current_bin <0)
      current_bin += h1->GetNbinsX();
    while (current_bin >= h1->GetNbinsX())
      current_bin -= h1->GetNbinsX();
    
    sum1 += h1->GetBinContent(current_bin+1);
    sum2 ++;
  }

  if (sum2>0){
    return sum1/sum2;
  }else{
    return 0;
  }
}

void plot_wROI(Int_t run=3493){
  Int_t plane = 0;
  //Int_t run = 3493;
  //  Int_t run = 5366;
  Int_t wire_bin_min; 
  Int_t wire_bin_max; 
  Int_t time_bin_min; 
  Int_t time_bin_max; 
  if (run == 3493){
    wire_bin_min = 1170; 
    wire_bin_max = 1300;
    // wire_bin_min = 3850-2400; 
    // wire_bin_max = 3930-2400;
    time_bin_min = 6200;
    time_bin_max = 7800;
  }else if (run==5366){
    wire_bin_min = 1220; 
    wire_bin_max = 1350;
    time_bin_min = 6300;
    time_bin_max = 7500;
  }

  Double_t factor1 = 6;
  time_bin_min /= factor1;
  time_bin_max /= factor1;

  TFile *file = new TFile(Form("%d_2D_all.root",run));
  TFile *file2 = new TFile(Form("%d_2D_lf_all_25.root",run));


  TH2F *hdecon;
  TH2F *hdecon1;
  TH2F *hdecon_th;
  if (plane == 0){
    hdecon = (TH2F*)file->Get("hdecon_u");
    hdecon_th = (TH2F*)file2->Get("hdecon_u");
  }else if (plane == 1){
    hdecon = (TH2F*)file->Get("hdecon_v");
    hdecon_th = (TH2F*)file2->Get("hdecon_v");
  }else if (plane == 2){
    hdecon = (TH2F*)file->Get("hdecon_w");
    hdecon_th = (TH2F*)file2->Get("hdecon_w");
  }
  hdecon->Rebin2D(1,factor1);
  hdecon_th->Rebin2D(1,factor1);

  const Int_t nwire = wire_bin_max - wire_bin_min + 1;
  Int_t ntime = time_bin_max - time_bin_min + 1;
  hdecon1 = new TH2F("hdecon1","hdecon1",nwire,0,nwire,ntime,0,ntime);

  TH1F *hwave1 = new TH1F("hwave1","hwave1",ntime,0,ntime);
  TH1F *hwave2 = new TH1F("hwave2","hwave2",ntime,0,ntime);
  Int_t special_channel = 51;
  Double_t factor = 3.5;
  Double_t factor2 = 0.7;
  Double_t max_th = 10000; // maximum threshold
  
  Double_t rms[nwire];
  Int_t ntime1 = hdecon_th->GetNbinsY();
  TH1F *htemp = new TH1F("htemp","htemp",ntime1,0,ntime1);
  for (Int_t i=wire_bin_min; i!=wire_bin_max+1;i++){
    for (Int_t j=0;j!=ntime1;j++){
      htemp->SetBinContent(j+1,hdecon_th->GetBinContent(i+1,j+1));
    }
    rms[i] = cal_rms1(htemp) * factor * factor2;
    if (rms[i] > max_th * factor2) rms[i] = max_th * factor2;
  }


  for (Int_t i=time_bin_min;i!=time_bin_max+1;i++){
    hwave1->SetBinContent(i-time_bin_min+1,hdecon->GetBinContent(wire_bin_min + special_channel+1,i+1));
  }
  
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
      if (s_end >= time_bin_min || s_begin <= time_bin_max){
	// get start and end
	Double_t content_begin = hdecon->GetBinContent(s_channel+1,s_begin+1);
	Double_t content_end = hdecon->GetBinContent(s_channel+1,s_end+1);
	
	for (Int_t j=s_begin; j!=s_end+1; j++){
	  Double_t content_current =  hdecon->GetBinContent(s_channel+1,j+1);
	  Double_t content = content_current - ((content_end - content_begin)*(j*1.-s_begin)/(s_end-s_begin*1.) + content_begin);
	  

	  if (j >= time_bin_min && j<= time_bin_max){
	    //	    if (content > rms[s_channel])
	      hdecon1->SetBinContent(s_channel-wire_bin_min+1,j-time_bin_min+1,
				     content);
	    if (s_channel == wire_bin_min + special_channel){
	      hwave2->SetBinContent(j-time_bin_min+1,content);
	      // std::cout << s_begin-time_bin_min << " " << s_end-time_bin_min << " " << content_current << " " << content << endl;
	    }
	  }
	  
	  


	}
      }
    }
  }
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  hdecon1->Draw("COLZ");
  c1->cd(2);
  hwave1->Draw();
  hwave2->SetLineColor(2);
  hwave2->Draw("same");
}
