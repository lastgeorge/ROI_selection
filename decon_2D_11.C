#include <vector>

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


void decon_2D_11(){
  Int_t plane = 1; // U plane
  TFile *file = new TFile("2D_display_3493_821_0.root");
  //TFile *file = new TFile("2D_display_5366_0_0.root");
  TH2I *hraw;
  // Int_t wire_bin_min = 1170; 
  // Int_t wire_bin_max = 1300;
  
  Int_t wire_bin_min = 3850; 
  Int_t wire_bin_max = 3930;
  Int_t time_bin_min = 6200;
  Int_t time_bin_max = 7800;
  

  // Int_t wire_bin_min = 1220; 
  // Int_t wire_bin_max = 1350;
  
  // Int_t time_bin_min = 6300;
  // Int_t time_bin_max = 7500;

  
  Int_t nwire = wire_bin_max - wire_bin_min + 1;
  Int_t ntime = time_bin_max - time_bin_min + 1;

  if (plane == 0){
    hraw = (TH2I*)file->Get("hu_raw");
  }else if (plane == 1){
    hraw = (TH2I*)file->Get("hv_raw");
    wire_bin_min -= 2400;
    wire_bin_max -= 2400;
  }else if (plane == 2){
    hraw = (TH2I*)file->Get("hw_raw");
    wire_bin_min -= 4800;
    wire_bin_max -= 4800;
  }
  
  // fill in the new content
  TH2F *hregion = new TH2F("hregion","hregion",nwire, wire_bin_min, wire_bin_max+1, ntime, time_bin_min, time_bin_max+1);
  TH2F *hdecon = new TH2F("hdecon","hdecon",nwire, wire_bin_min, wire_bin_max+1, ntime, time_bin_min, time_bin_max+1);
  TH2F *hregiona = new TH2F("hregiona","hregiona",nwire, wire_bin_min, wire_bin_max+1, ntime, time_bin_min, time_bin_max+1);
  for (Int_t i=wire_bin_min; i!= wire_bin_max+1; i++){
    for (Int_t j=time_bin_min; j!=time_bin_max+1; j++){
      hregion->SetBinContent(i-wire_bin_min+1,j-time_bin_min+1,hraw->GetBinContent(i+1,j+1));
      hregiona->SetBinContent(i-wire_bin_min+1,j-time_bin_min+1,fabs(hraw->GetBinContent(i+1,j+1)));
    }
  }
  //hregion->Draw("COLZ");
  
  // read in the field responses
  // which plane
#include "data_70_ROI_11.txt"
  TGraph **gf_2D_g = new TGraph*[11];

  // for (Int_t i=0;i!=5000;i++){
  //   u_2D_g_4_y[i] = 0;
  //   u_2D_g_5_y[i] = 0;
  //   u_2D_g_6_y[i] = 0;
  //   u_2D_g_7_y[i] = 0;
  //   u_2D_g_8_y[i] = 0;
  //   u_2D_g_9_y[i] = 0;
  //   u_2D_g_10_y[i] = 0;
  // }

  if (plane == 0){
    gf_2D_g[0] = new TGraph(5000,u_2D_g_0_x,u_2D_g_0_y);
    gf_2D_g[1] = new TGraph(5000,u_2D_g_1_x,u_2D_g_1_y);
    gf_2D_g[2] = new TGraph(5000,u_2D_g_2_x,u_2D_g_2_y);
    gf_2D_g[3] = new TGraph(5000,u_2D_g_3_x,u_2D_g_3_y);
    gf_2D_g[4] = new TGraph(5000,u_2D_g_4_x,u_2D_g_4_y);
    gf_2D_g[5] = new TGraph(5000,u_2D_g_5_x,u_2D_g_5_y);
    gf_2D_g[6] = new TGraph(5000,u_2D_g_6_x,u_2D_g_6_y);
    gf_2D_g[7] = new TGraph(5000,u_2D_g_7_x,u_2D_g_7_y);
    gf_2D_g[8] = new TGraph(5000,u_2D_g_8_x,u_2D_g_8_y);
    gf_2D_g[9] = new TGraph(5000,u_2D_g_9_x,u_2D_g_9_y);
    gf_2D_g[10] = new TGraph(5000,u_2D_g_10_x,u_2D_g_10_y);
  }else if (plane == 1){
    gf_2D_g[0] = new TGraph(5000,v_2D_g_0_x,v_2D_g_0_y);
    gf_2D_g[1] = new TGraph(5000,v_2D_g_1_x,v_2D_g_1_y);
    gf_2D_g[2] = new TGraph(5000,v_2D_g_2_x,v_2D_g_2_y);
    gf_2D_g[3] = new TGraph(5000,v_2D_g_3_x,v_2D_g_3_y);
    gf_2D_g[4] = new TGraph(5000,v_2D_g_4_x,v_2D_g_4_y);
    gf_2D_g[5] = new TGraph(5000,v_2D_g_5_x,v_2D_g_5_y);
    gf_2D_g[6] = new TGraph(5000,v_2D_g_6_x,v_2D_g_6_y);
    gf_2D_g[7] = new TGraph(5000,v_2D_g_7_x,v_2D_g_7_y);
    gf_2D_g[8] = new TGraph(5000,v_2D_g_8_x,v_2D_g_8_y);
    gf_2D_g[9] = new TGraph(5000,v_2D_g_9_x,v_2D_g_9_y);
    gf_2D_g[10] = new TGraph(5000,v_2D_g_10_x,v_2D_g_10_y);
  }
  // prepare for resulse functions ...
  const Int_t nticks = ntime;
  const Int_t nchannels = nwire;

  double rho_res[21][nticks], phi_res[21][nticks];

  TF1 *filter_wire = new TF1("filter_wire","exp(-0.5*pow(x/[0],2))");
  double par4[1] = {1.0/sqrt(3.1415926)*1.4};
  filter_wire->SetParameters(par4);
  TF1 *filter_time = new TF1("filter_u","(x>0.0)*exp(-0.5*pow(x/[0],[1]))");
  double par[2]={1.43555e+01/200.*2.,4.95096e+00};
  filter_time->SetParameters(par);
  //hregion->SetBinContent(i-wire_bin_min+1,j-time_bin_min+1,hraw->GetBinContent(i+1,j+1));
  //TF1 *filter_low = new TF1("filter_low","1-exp(-pow(x/0.0045,2))");
  TF1 *filter_low = new TF1("filter_low","1-exp(-pow(x/0.0025,2))");

  TH1F *hvr = new TH1F("hvr1","hvr1",nticks,0,nticks); // half us tick
  for (int j=0;j!=21;j++){
    TGraph *gtemp;
    if (j==0 || j==20){
      gtemp = gf_2D_g[10];
    }else if (j==1 || j==19){
      gtemp = gf_2D_g[9];
    }else if (j==2 || j==18){
      gtemp = gf_2D_g[8];
    }else if (j==3 || j==17){
      gtemp = gf_2D_g[7];
    }else if (j==4 || j==16){
      gtemp = gf_2D_g[6];
    }else if (j==5 || j==15){
      gtemp = gf_2D_g[5];
    }else if (j==6 || j==14){
      gtemp = gf_2D_g[4];
    }else if (j==7 || j==13){
      gtemp = gf_2D_g[3];
    }else if (j==8 || j==12){
      gtemp = gf_2D_g[2];
    }else if (j==9 || j==11){
      gtemp = gf_2D_g[1];
    }else if (j==10){
      gtemp = gf_2D_g[0];
    }
    hvr->Reset();
    for (int i=0; i!=nticks; i++){  
      double time = hvr->GetBinCenter(i+1)/2.- 90 ;
      //*** scale factors for 70kV ***//
      double x = time; // 0.113 is the calibration number ... 
      if (x > -84 && x  < 15.8){
	hvr->SetBinContent(i+1,gtemp->Eval(x)); //70kV
      }else{
	hvr->SetBinContent(i+1,0);
      }
    }
    TH1 *hmr_v = hvr->FFT(0,"MAG");
    TH1 *hpr_v = hvr->FFT(0,"PH");
    
    for (Int_t i=0;i!=nticks;i++){
      rho_res[j][i] = hmr_v->GetBinContent(i+1);
      phi_res[j][i] = hpr_v->GetBinContent(i+1);
    }

    // std::cout << "abc " << j << std::endl;
    delete hmr_v;
    delete hpr_v;
  }

  // do the 2D deconvolution
  std::vector<std::vector<double>> rho_v, phi_v,result_re,result_im;
  for (int i=0;i!=nchannels;i++){
    std::vector<double> temp,temp1,temp2,temp3;
    temp.resize(nticks,0);
    rho_v.push_back(temp);
    temp1.resize(nticks,0);
    phi_v.push_back(temp1);
    temp2.resize(nticks,0);
    result_re.push_back(temp2);
    temp3.resize(nticks,0);
    result_im.push_back(temp3);
  }
  int tbin_save[nchannels];
  
  for (Int_t ind = 0; ind!= nchannels; ind ++){
    TH1F *htemp = new TH1F("htemp","htemp",nticks,0,nticks);
    for (Int_t i=0;i!=nticks;i++){
      htemp->SetBinContent(i+1,hregion->GetBinContent(ind+1,i+1));
    }
    TH1 *hm = htemp->FFT(0,"MAG");
    TH1 *hp = htemp->FFT(0,"PH");

    for (Int_t j=0;j!=nticks;j++){
      rho_v[ind][j] = hm->GetBinContent(j+1);
      phi_v[ind][j] = hp->GetBinContent(j+1);
    }
    tbin_save[ind] = 0;
    delete hm;
    delete hp;
    delete htemp;
  }

  double resp_re[nchannels], resp_im[nchannels];
  double value_re[9600];
  double value_im[9600];

  for (Int_t i=0;i!=nticks;i++){
    Double_t freq;
    if (i < nticks/2.){
      freq = i/(1.*nticks)*2.;
    }else{
      freq = (nticks - i)/(1.*nticks)*2.;
    }
    
    for (Int_t j=0;j!=nchannels;j++){
      value_re[j] = rho_v[j][i]*cos(phi_v[j][i]);
      value_im[j] = rho_v[j][i]*sin(phi_v[j][i]);
      if (j<21){
  	resp_re[j] = rho_res[j][i]*cos(phi_res[j][i]);
  	resp_im[j] = rho_res[j][i]*sin(phi_res[j][i]);
      }else{
  	resp_re[j] = 0.;
  	resp_im[j] = 0.;
      }
    }

    Int_t m=nchannels;
    
    TVirtualFFT *ifft = TVirtualFFT::FFT(1,&m,"C2CFORWARD M K");
    ifft->SetPointsComplex(value_re,value_im);
    ifft->Transform();
    Double_t temp_re[nchannels],temp_im[nchannels];
    ifft->GetPointsComplex(temp_re,temp_im);
    
    ifft->SetPointsComplex(resp_re,resp_im);
    ifft->Transform();
    Double_t temp1_re[nchannels],temp1_im[nchannels];
    ifft->GetPointsComplex(temp1_re,temp1_im);
    
    Double_t temp2_re[nchannels],temp2_im[nchannels];
    for (Int_t j=0;j!=nchannels;j++){
      Double_t freq_wire;
      if (j < nchannels/2.){
	freq_wire = j/(1.*nchannels)*2.;
      }else{
	freq_wire = (nchannels -j)/(1.*nchannels)*2.;
      }
      if (temp1_im[j]*temp1_im[j]+temp1_re[j]*temp1_re[j]>0){
  	temp2_re[j] = (temp_re[j]*temp1_re[j]+temp_im[j]*temp1_im[j])/m/
  	  (temp1_im[j]*temp1_im[j]+temp1_re[j]*temp1_re[j]) * filter_wire->Eval(freq_wire);
  	temp2_im[j] = (temp_im[j]*temp1_re[j]-temp_re[j]*temp1_im[j])
  	  /m/(temp1_im[j]*temp1_im[j]+temp1_re[j]*temp1_re[j]) * filter_wire->Eval(freq_wire);
      }else{
  	temp2_re[j] = 0;
  	temp2_im[j] = 0;
      }
    }
    
    TVirtualFFT *ifft3 = TVirtualFFT::FFT(1,&m,"C2CBACKWARD M K");
    ifft3->SetPointsComplex(temp2_re,temp2_im);
    ifft3->Transform();
    Double_t temp3_re[nchannels],temp3_im[nchannels];
    ifft3->GetPointsComplex(temp3_re,temp3_im);


    for (Int_t j=0;j!=nchannels;j++){
      Int_t shift = j - 10;
      if (shift <0) shift += nchannels;
      result_re[j][i] = temp3_re[shift]/nticks;
      result_im[j][i] = temp3_im[shift]/nticks;
    }

    delete ifft;
    delete ifft3;
  }

  for (Int_t chid=0;chid!=nchannels;chid++){
    int n = nticks;
    TVirtualFFT *ifft2 = TVirtualFFT::FFT(1,&n,"C2R M K");
    double temp_re[nticks],temp_im[nticks];
    for (int j=0;j!=nticks;j++){
      Double_t freq;
      if (j < nticks/2.){
	freq = j/(1.*nticks)*2.;
      }else{
	freq = (nticks - j)/(1.*nticks)*2.;
      }
      temp_re[j] = result_re[chid][j]*filter_time->Eval(freq) * filter_low->Eval(freq);
      temp_im[j] = result_im[chid][j]*filter_time->Eval(freq) * filter_low->Eval(freq);
    }
    ifft2->SetPointsComplex(temp_re,temp_im);
    ifft2->Transform();
    TH1 *fb = 0;
    fb = TH1::TransformHisto(ifft2,fb,"Re");
    
    TH1F *htemp = new TH1F("htemp","htemp",nticks,0,nticks);
    for (int i=0;i!=nticks;i++){
      htemp->SetBinContent(i+1,fb->GetBinContent(i+1)/( 14.*1.1*4096./2000.));
    }
    delete ifft2;
    delete fb;
    
    // correct baseline 
    restore_baseline(htemp);
    
    // put results back into the 2-D histogram
    for (int j=0;j!=nticks;j++){
      int bin = j+180;
      if (bin >= nticks) bin -= nticks;
      hdecon->SetBinContent(chid+1,bin+1,htemp->GetBinContent(j+1));
    }
    delete htemp;
   }

  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  //c1_1->SetLogz(1);
  hregiona->Draw("COLZ");
  c1->cd(2);
  //c1_2->SetLogz(1);
  hdecon->Draw("COLZ");

  hregion->SetTitle("Raw");
  hdecon->SetTitle("2D Decon");

    // save the images
}
