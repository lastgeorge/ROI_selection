Double_t PoissonReal(const Double_t *k,  const Double_t *lambda) {
  return TMath::Exp(k[0]*TMath::Log(lambda[0])-lambda[0]) / TMath::Gamma(k[0]+1.);
}

void do_toy_mc(){
  TFile *file1 = new TFile("2D_display_3493_821_0.root");
  TH2F *hu_raw = (TH2F*)file1->Get("hu_raw");
  TH2F *hv_raw = (TH2F*)file1->Get("hv_raw");
  TH2F *hw_raw = (TH2F*)file1->Get("hw_raw");

  hu_raw->Reset();
  hv_raw->Reset();
  hw_raw->Reset();

  Double_t MaxPoissArg = 100.; 
  
  TF1 *MyPoisson = new TF1("MyPoisson", PoissonReal,0.,MaxPoissArg,1);
  //**Setting lambda**//
  Double_t params[1] = {3.3708}; //2us
  MyPoisson->SetParameters(params);

  TF1 *f1 = new TF1("f","[0]+([1]+[2]*x*9592./2.)*exp(-[3]*pow(x*9592./2.,[4]))"); // x in MHz
  Double_t fitpar[5]={4.27132e+01,6.22750e+02,-2.53535e-01,8.07578e-05,1.35510e+00}; //2us (Final Fit)
  f1->SetParameters(fitpar);
  

  Double_t value_re[10000], value_im[10000];
  Int_t  n  = hu_raw->GetNbinsY(); // window size 

  for (Int_t k=0;k!=hu_raw->GetNbinsX();k++){
    if (k%100==0) cout << "U " << k << endl;
     for (Int_t i=0;i!=n;i++){
      Double_t freq;
      if (i<n/2.){
	freq = (i)*2./n; // assume 2 MHz digitization 
      }else{
	freq = (n-i)*2./n;
      }
      
      Double_t rho = f1->Eval(freq) * sqrt(n/9592.) * MyPoisson->GetRandom()/params[0];
      if (i==0) rho = 0;
      Double_t phi = gRandom->Uniform(0,2*3.1415926);
      value_re[i] = rho*cos(phi)/n;
      value_im[i] = rho*sin(phi)/n;
    }
    
    TVirtualFFT *ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
    ifft->SetPointsComplex(value_re,value_im);
    ifft->Transform();
    TH1 *fb = 0;
    fb = TH1::TransformHisto(ifft,fb,"Re");
    // fb->Draw();
    
    fb->Scale(2.1/2.0);
    
    for (Int_t i=0;i!=n;i++){
      hu_raw->SetBinContent(k+1,i+1,fb->GetBinContent(i+1));
    }
    delete fb;
  }


   for (Int_t k=0;k!=hv_raw->GetNbinsX();k++){
    if (k%100==0) cout << "V " << k << endl;
     for (Int_t i=0;i!=n;i++){
      Double_t freq;
      if (i<n/2.){
	freq = (i)*2./n; // assume 2 MHz digitization 
      }else{
	freq = (n-i)*2./n;
      }
      
      Double_t rho = f1->Eval(freq) * sqrt(n/9592.) * MyPoisson->GetRandom()/params[0];
      if (i==0) rho = 0;
      Double_t phi = gRandom->Uniform(0,2*3.1415926);
      value_re[i] = rho*cos(phi)/n;
      value_im[i] = rho*sin(phi)/n;
    }
    
    TVirtualFFT *ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
    ifft->SetPointsComplex(value_re,value_im);
    ifft->Transform();
    TH1 *fb = 0;
    fb = TH1::TransformHisto(ifft,fb,"Re");
    // fb->Draw();
    
    fb->Scale(2.1/2.0);
    
    for (Int_t i=0;i!=n;i++){
      hv_raw->SetBinContent(k+1,i+1,fb->GetBinContent(i+1));
    }
    delete fb;
  }


   for (Int_t k=0;k!=hw_raw->GetNbinsX();k++){
    if (k%100==0) cout << "W " << k << endl;
     for (Int_t i=0;i!=n;i++){
      Double_t freq;
      if (i<n/2.){
	freq = (i)*2./n; // assume 2 MHz digitization 
      }else{
	freq = (n-i)*2./n;
      }
      
      Double_t rho = f1->Eval(freq) * sqrt(n/9592.) * MyPoisson->GetRandom()/params[0];
      if (i==0) rho = 0;
      Double_t phi = gRandom->Uniform(0,2*3.1415926);
      value_re[i] = rho*cos(phi)/n;
      value_im[i] = rho*sin(phi)/n;
    }
    
    TVirtualFFT *ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
    ifft->SetPointsComplex(value_re,value_im);
    ifft->Transform();
    TH1 *fb = 0;
    fb = TH1::TransformHisto(ifft,fb,"Re");
    // fb->Draw();
    
    fb->Scale(1.5/2.0);
    
    for (Int_t i=0;i!=n;i++){
      hw_raw->SetBinContent(k+1,i+1,fb->GetBinContent(i+1));
    }
    delete fb;
  }



  TFile *file2 = new TFile("toy_mc.root","RECREATE");
  hu_raw->SetDirectory(file2);
  hv_raw->SetDirectory(file2);
  hw_raw->SetDirectory(file2);
  file2->Write();
  file2->Close();

  //cout << hu_raw->GetNbinsY() << endl;

}
