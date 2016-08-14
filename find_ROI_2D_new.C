#include <vector>
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include <iostream>

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


Int_t find_ROI_end(TH1F *h1, Int_t bin, Double_t th=0){
  Int_t end = bin;
  Double_t content = h1->GetBinContent(end+1);
  while(content>th){
    end ++;
    if (end >=h1->GetNbinsX()){
      content = h1->GetBinContent(end - h1->GetNbinsX()+1);
    }else{
      content = h1->GetBinContent(end+1);
    }
    if (end == h1->GetNbinsX()) break;
  }

  while(local_ave(h1,end+1,1) < local_ave(h1,end,1)){
    end++;
    if (end == h1->GetNbinsX()) break;
  } 
  return end;
}

Int_t find_ROI_begin(TH1F *h1, Int_t bin, Double_t th=0){
  // find the first one before bin and is below threshold ... 
  Int_t begin = bin;
  Double_t content = h1->GetBinContent(begin+1);
  while(content > th){
    begin --;
    if (begin <0){
      content = h1->GetBinContent(begin+h1->GetNbinsX()+1);
    }else{
      content = h1->GetBinContent(begin+1);
    }
    if (begin == 0) break;
  }
  
  // calculate the local average
  // keep going and find the minimum
  while( local_ave(h1,begin-1,1) < local_ave(h1,begin,1)){
    begin --;
    if (begin == 0) break;
  }
  
  return begin;
}





void find_ROI_2D_new(Int_t run = 3455){
  TFile *file = new TFile(Form("%d_2D_lf_all_25.root",run));
  TH2F *hdecon_u = (TH2F*)file->Get("hdecon_u");
  TH2F *hdecon_v = (TH2F*)file->Get("hdecon_v");
  TH2F *hdecon_w = (TH2F*)file->Get("hdecon_w");
  
  //rebin
  hdecon_u->Rebin2D(1,6);
  hdecon_v->Rebin2D(1,6);
  hdecon_w->Rebin2D(1,6);
  
  //save ROIs
  TFile *file1 = new TFile(Form("ROI_2D_25_%d.root",run),"RECREATE");
  TTree *T = new TTree("T","T");
  T->SetDirectory(file1);
  Int_t s_begin,s_end,s_channel,s_plane;
  T->Branch("begin",&s_begin,"begin/I");
  T->Branch("end",&s_end,"end/I");
  T->Branch("channel",&s_channel,"channel/I");
  T->Branch("plane",&s_plane,"plane/I");
  Int_t ntime = hdecon_u->GetNbinsY();
  TH1F *htemp = new TH1F("htemp","htemp",ntime,0,ntime);
  Int_t nchannels = hdecon_u->GetNbinsX();
  

  // some example thresholds
  Double_t factor = 3.5; // regular threshold
  Double_t max_th = 10000; // maximum threshold
  Double_t factor1 = 0.7; //special threshold
  Int_t short_length = 3; // short length
  Double_t fixed_threshold = 4000;
  
  s_plane = 0;
  for (Int_t i=0;i!=nchannels;i++){
    //for (Int_t i=1282;i!=1283;i++){
    s_channel = i;
    if (i%100 ==0) std::cout << i << std::endl;

    for (Int_t j=0;j!=ntime;j++){
      htemp->SetBinContent(j+1,hdecon_u->GetBinContent(i+1,j+1));
    }

    //htemp->Draw();
    // std::vector<std::pair <Int_t,Int_t> > ROIs;
    std::vector<std::pair <Int_t,Int_t> > ROIs_1;
    //    std::vector<Int_t> max_bins;
    std::vector<Int_t> max_bins_1;

    Double_t th = cal_rms1(htemp) * factor;
    if (th > max_th) th = max_th;

    // first round ROIs
    for (Int_t j=1; j<ntime-1;j++){
      Double_t content = htemp->GetBinContent(j+1);
      Double_t prev_content = htemp->GetBinContent(j);
      Double_t next_content = htemp->GetBinContent(j+2);
      Int_t flag_ROI = 0;
      Int_t begin;
      Int_t end;
      Int_t max_bin;

      if (content > th){
      	begin = find_ROI_begin(htemp,j, th*factor1) ;
	end = find_ROI_end(htemp,j, th*factor1) ;
	max_bin = begin;
	for (Int_t k=begin;k<=end;k++){
	  if (htemp->GetBinContent(k+1) > htemp->GetBinContent(max_bin+1)){
	    max_bin = k;
	  }
	}
	flag_ROI = 1;
      }else{
	if (content > prev_content && content > next_content){
	  begin = find_ROI_begin(htemp,j, prev_content);
	  end = find_ROI_end(htemp,j, next_content );
	  max_bin = begin;
	  for (Int_t k=begin;k<=end;k++){
	    if (htemp->GetBinContent(k+1) > htemp->GetBinContent(max_bin+1)){
	      max_bin = k;
	    }
	  }
	  
	  // if the difference is big enough
	  if (htemp->GetBinContent(max_bin+1) - htemp->GetBinContent(begin+1) +
	      htemp->GetBinContent(max_bin+1) - htemp->GetBinContent(end+1) > th*2*factor1){
	    flag_ROI = 1;
	  }

	  //	  cout << max_bin << " " << begin << " " << end << " " << htemp->GetBinContent(max_bin+1) << " " << htemp->GetBinContent(begin+1) << " " << htemp->GetBinContent(end+1) << endl;
	  // if the whole thing is quite narrow
	  
	  Int_t temp_begin = max_bin-short_length;
	  if (temp_begin < begin) temp_begin = begin;
	  Int_t temp_end = max_bin + short_length;
	  if (temp_end > end) temp_end = end;
	  if ((htemp->GetBinContent(max_bin+1) - htemp->GetBinContent(temp_begin+1) > th *factor1 &&
	       htemp->GetBinContent(max_bin+1) - htemp->GetBinContent(temp_end+1) > th*factor1) ){
	    flag_ROI = 1;
	  }
	  
	}
      }

      if (flag_ROI == 1){
	if (ROIs_1.size()>0){
	  if (begin <= ROIs_1.back().second){
	    ROIs_1.back().second = end;
	    if (htemp->GetBinContent(max_bin+1) > htemp->GetBinContent(max_bins_1.back()+1))
	      max_bins_1.back() = max_bin;
	  }else{
	    ROIs_1.push_back(std::make_pair(begin,end));
	    max_bins_1.push_back(max_bin);
	  }
	}else{
	  ROIs_1.push_back(std::make_pair(begin,end));
	  max_bins_1.push_back(max_bin);
	}
	
	if (end <htemp->GetNbinsX()){
	  j = end;
	}else{
	  j=htemp->GetNbinsX();
	}
      }
    }
       
    // // Now go through all the gaps between ROIs, consider if need to merge them ... this is to take care of broad peak ... 
    // // find the minimum between the two ROIs, draw a line, and evaluate the middle part, make sure everybody pass threshold ... 
    if (ROIs_1.size()==1){
    }else if (ROIs_1.size()>1){
      Int_t flag_repeat = 1;
      //  cout << "Xin1: " << ROIs_1.size() << endl;;
      while(flag_repeat){
	flag_repeat = 0;
	for (Int_t k=0;k<ROIs_1.size()-1;k++){
	  Int_t begin = ROIs_1.at(k).first;
	  Int_t end = ROIs_1.at(k+1).second;
	  
	  Double_t begin_content = htemp->GetBinContent(begin+1);
	  Double_t end_content = htemp->GetBinContent(end+1);
	  
	  Int_t begin_1 = ROIs_1.at(k).second;
	  Int_t end_1 = ROIs_1.at(k+1).first;	
	  
	  Int_t flag_merge = 1;
	  //Double_t sum1 = 0, sum2 = 0;
	  for (Int_t j=begin_1; j<=end_1;j++){
	    Double_t current_content = htemp->GetBinContent(j+1);
	    Double_t content = current_content - ((end_content - begin_content)*(j*1.-begin)/(end-begin*1.) + begin_content);
	    
	    if (content < th*factor1){
	      flag_merge = 0;
	      break;
	    }
	    // sum1 += content;
	    // sum2 ++;
	    // cout << j << " " << content << endl;
	  }
	  // if (sum2 >0){
	  //   if (sum1/sum2 < th*factor1) flag_merge = 0;
	  // }
	  
	  if (flag_merge == 1){
	    ROIs_1.at(k).second = ROIs_1.at(k+1).second;
	    ROIs_1.erase(ROIs_1.begin()+k+1);
	    flag_repeat = 1;
	    break;
	  }
	}
	//	cout << "Xin2: " << ROIs_1.size() << endl;

      }
    }
    

    // sort out the ROI and save ... 
    for (auto it = ROIs_1.begin();it!=ROIs_1.end();it++){
      s_begin = it->first;
      s_end = it->second;
      //      cout << s_begin << " " << s_end << " " << th << " " << th *factor1 << endl;
      T->Fill();
    }
    
  }


  htemp->Draw();
 

  
  file1->Write();
  file1->Close();

}
