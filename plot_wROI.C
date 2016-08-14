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
  TH2F *hdecon;
  TH2F *hdecon1;
  if (plane == 0){
    hdecon = (TH2F*)file->Get("hdecon_u");
  }else if (plane == 1){
    hdecon = (TH2F*)file->Get("hdecon_v");
  }else if (plane == 2){
    hdecon = (TH2F*)file->Get("hdecon_w");
  }
  hdecon->Rebin2D(1,factor1);

  Int_t nwire = wire_bin_max - wire_bin_min + 1;
  Int_t ntime = time_bin_max - time_bin_min + 1;
  hdecon1 = new TH2F("hdecon1","hdecon1",nwire,0,nwire,ntime,0,ntime);

  TH1F *hwave1 = new TH1F("hwave1","hwave1",ntime,0,ntime);
  TH1F *hwave2 = new TH1F("hwave2","hwave2",ntime,0,ntime);
  Int_t special_channel = 101;


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
	    if (content > 3000)
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
