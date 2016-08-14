void plot_rebin(){
  ifstream infile1("result_r1.dat");
  ifstream infile2("result_r2.dat");
  ifstream infile3("result_r4.dat");
  ifstream infile4("result_r6.dat");
  Double_t x,y;
  Double_t temp;
  TGraph *graw = new TGraph();
  TGraph *g1D = new TGraph();
  TGraph *g2D_lf_25 = new TGraph();
  TGraph *g2D_lf_45 = new TGraph();
  for (Int_t i=0;i!=8256;i++){
    x = i+1;

    infile1 >> temp >> y;
    // y *= 198;    
    graw->SetPoint(i,x,y);
    
    infile2 >> temp >> y;
    //y*=5;
    g1D->SetPoint(i,x,y);
    
    infile3 >> temp >> y;
    //  y*=2.5;
    g2D_lf_25->SetPoint(i,x,y);

    infile4 >> temp >> y;
    //  y*=2.5;
    g2D_lf_45->SetPoint(i,x,y);
  }

  
  graw->Draw("A*");
  graw->Fit("pol0","","",3000,4000);
  graw->Fit("pol0","","",5000,8000);

  g1D->Draw("*same");
  g1D->Fit("pol0","","",3000,4000);
  g1D->Fit("pol0","","",5000,8000);
  g2D_lf_25->Draw("*same");
  g2D_lf_25->Fit("pol0","","",3000,4000);
  g2D_lf_25->Fit("pol0","","",5000,8000);
  g2D_lf_45->Draw("*same");
  g2D_lf_45->Fit("pol0","","",3000,4000);
  g2D_lf_45->Fit("pol0","","",5000,8000);

  graw->SetMarkerStyle(20);
  g1D->SetMarkerStyle(20);
  g2D_lf_25->SetMarkerStyle(20);
  g2D_lf_45->SetMarkerStyle(20);
  
  graw->SetMarkerColor(1);
  g1D->SetMarkerColor(2);
  g2D_lf_25->SetMarkerColor(4);
  g2D_lf_45->SetMarkerColor(6);
}
