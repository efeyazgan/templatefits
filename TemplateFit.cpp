#include <math.h>
#include <iostream>
#include <algorithm>
#include <string.h>
#include <sstream>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include "TROOT.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TFile.h"
#include "RooPolyVar.h"
#include "TStyle.h"
#include "RooDataHist.h"
#include <RooDataSet.h>
#include <RooHistPdf.h>
#include <RooGaussian.h>
#include <RooRealVar.h>
#include <RooMinuit.h>
#include "RooMinimizer.h"
#include <RooFitResult.h>
#include <RooGlobalFunc.h>
#include <RooPlot.h>
#include "RooTFnBinding.h"
#include <TMultiGraph.h>
#include "TF2.h"
#include "RooGenericPdf.h"
#include "RooCategory.h"
#include "RooEfficiency.h"
#include "RooPlot.h"
#include "TLegend.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "TROOT.h"
#include "RooGlobalFunc.h" 
#include <stdio.h>
#include "TStyle.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include "TFile.h"
#include "TLegend.h"
#include "RooRealVar.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooCategory.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TVector.h"
#include "TCanvas.h"
#include "TProfile.h"
#include <math.h>


float makeArray( string  INFILE_NAME, int Ev, bool uncor){

  stringstream name;
  name<<"/data/weights/"<<INFILE_NAME<<"/weights-"<<Ev<<".out";
  stringstream name2;
  int Ev2 = Ev+50;
  name2<<"/data/weights/"<<INFILE_NAME<<"/weights-"<<Ev2<<".out";
  //Open input file                                                                                                                                      
  ifstream myInfile(name.str());
  ifstream myInfile2(name2.str());
  string event_number;
  string likelihood;
  string likelihood_unc;
  float L, L_unc;
  float Like = 0;
  float Like_unc = 0;

  if ( myInfile.is_open() ) {
    while ( myInfile>>event_number){
      myInfile>>likelihood; myInfile>>likelihood_unc;
      if(uncor){L = log(atof(likelihood.c_str()));}
      else{L = log(atof(likelihood.c_str()));}
      L_unc = log(atof(likelihood_unc.c_str()));
      Like += L;
      Like_unc += L_unc*L_unc;
    }
  }
  myInfile.close();

  Like_unc = sqrt(Like_unc);


  return Like;

}

float SampleLikelihoodH(std::vector<float> *Cor, std::vector<float> *UnCor, float H)
{
  std::vector<float> *Event_L = new std::vector<float>();
  float Likelihood = 0.0;
  for(unsigned int i = 0; i < Cor->size(); i++){
    Event_L->push_back(-log(H * (*Cor)[i] + (1.0 - H) * (*UnCor)[i]));
    Likelihood += (-1)*log( H * (*Cor)[i] + (1.0 - H) * (*UnCor)[i]);    
  }
  return Likelihood;
}

void MakeGraph(float X[], float Y[], TString name)
{
  TGraph *GraphCase1 = new TGraph(11,X,Y);
  GraphCase1->SetTitle("Sample Likelihood");
  GraphCase1->GetXaxis()->SetTitle("H");
  GraphCase1->GetYaxis()->SetTitle("Likelihood");
  GraphCase1->SetMarkerStyle(21);
  GraphCase1->SetMarkerSize(1);
  TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,500);
  c1->Range(0,0,1,1);
  c1->cd();
  GraphCase1->Draw("AP");
  c1->Modified();
  c1->Update();
  c1->SaveAs( name);

}

std::vector<double> DoFit(TH1F *hist_data, TH1F* hist_MC_cor, TH1F* hist_MC_uncor){

  gSystem->Load("libRooFit") ;
  
  
  // The fit variable - lepton invariant mass                                                                                                                                                                
  RooRealVar* IndL_ = new RooRealVar("IndL","IndL",-3., 3.0, "");
  RooRealVar IndL = *IndL_;
  
  gROOT->cd();
  
  RooArgList X(IndL);
  ///////// convert Histograms into RooDataHists                                                                                                           
  RooDataHist* data = new RooDataHist("data","data",X, hist_data,2);
  RooDataHist* MC_cor = new RooDataHist("MC_cor","MC_cor", X, hist_MC_cor,2);
  RooDataHist* MC_uncor = new RooDataHist("MC_uncor","MC_uncor", X, hist_MC_uncor,2);
  
  //  RooRealVar* Cor_mean = new RooRealVar("Cor_mean","Cor_mean",0.238928);
  //RooRealVar* Cor_sigma = new RooRealVar("Cor_sigma","Cor_sigma",0.241596);
  //RooRealVar* UnCor_mean = new RooRealVar("UnCor_mean","UnCor_mean",0.193420);
  //RooRealVar* UnCor_sigma = new RooRealVar("UnCor_sigma","UnCor_sigma",0.249401);
  //RooGaussian* Cor = new RooGaussian("Cor","Cor",IndL,*Cor_mean, *Cor_sigma);
  //RooGaussian* UnCor = new RooGaussian("UnCor","UnCor",IndL,*UnCor_mean,*UnCor_sigma);
  RooHistPdf* Cor = new RooHistPdf("Cor","Cor", X, *MC_cor);
  RooHistPdf* UnCor = new RooHistPdf("UnCor","UnCor", X, *MC_uncor);
  
  // Now define some efficiency/yield variables                                                                                                                                                              
  RooRealVar* numSignal = new RooRealVar("numSignal","numSignal", 1, -1, 1);
  RooRealVar* numBackground = new RooRealVar("numBackground","numBackground", 0, -1, 1);
  
  //  RooArgList Sum(*Cor, *UnCor);
  RooArgList *Sum = new RooArgList();
  Sum->add(*Cor);
  Sum->add(*UnCor);
  RooArgList yieldsPass(*numSignal,*numBackground);
  //  RooRealVar* f = new RooRealVar("f","f",hist_data->Integral(),0.,20000.);
  RooRealVar* f = new RooRealVar("f","f",0.5,0.,1.);
  //    RooRealVar* fb = new RooRealVar("fb","fb",0.5,0.,1.);
  //RooArgList F(*f, *fb);
  //  RooAddPdf pdfPass("pdfPass","extended sum pdf", RooArgList( *Cor, *UnCor), RooArgList( *f, *fb));
  RooAddPdf pdfPass("pdfPass","extended sum pdf", *Sum, *f);

  // ********* Do the Actual Fit ********** //                                                                                                                                                               
  RooFitResult *fitResult = pdfPass.fitTo(*data, RooFit::Extended(1));
  TString cname = TString("fit");
  TCanvas* c = new TCanvas(cname,cname,500,500);
  RooPlot* frame1 = IndL.frame();
  frame1->SetMinimum(0);
  //  data->plotOn(frame1,RooFit::DataError(errorType));
  pdfPass.plotOn(frame1,RooFit::ProjWData(*data),RooFit::Components(*Cor),RooFit::LineColor(kRed));
  pdfPass.plotOn(frame1,RooFit::ProjWData(*data),RooFit::Components(*UnCor),RooFit::LineColor(kGreen));
  pdfPass.plotOn(frame1,RooFit::ProjWData(*data));
  // pdfPass.plotOn(frame1);
  pdfPass.paramOn(frame1,data);
  frame1->Draw("e0");
  c->SaveAs("fit.C");
  delete c;

  std::vector<double> results;
  //  double x = f->getVal()/(f->getVal() + fb->getVal());
  double x = f->getVal();
  results.push_back(x);
  //double Err_A = (f->getError()/f->getVal());
  //double Err_B = sqrt(f->getError()*f->getError() + fb->getError()*fb->getError())/(f->getVal()+fb->getVal());
  //double Err = x*sqrt(Err_A*Err_A + Err_B*Err_B);
  double Err = f->getError();
  results.push_back(Err);
  //  std::cout<<"numSignal= "<<numSignal->getVal()<<" numBackground= "<<numBackground->getVal()<<" f = numSignal/(numSignal + numBackground) = "<<numSignal->getVal()/(numSignal->getVal() + numBackground->getVal())<<std::endl;
  std::cout<<"x= "<<(results)[0]<<" , x_err= "<<(results)[1]<<std::endl;
  return results;
}

int main ( int argc, const char* argv[] ) {

  using namespace RooFit ;

  TH1F* GaussEvCor = new TH1F("Correlated","-2ln#lambda = 2lnL_{H=C} - 2lnL_{H=U}", 100, -3, 3);
  TH1F* GaussEvUnCor = new TH1F("UnCorrelated","-2ln#lambda = 2lnL_{H=C} - 2lnL_{H=U}", 100, -3, 3);
  
  for(int i=0; i < 200; i++){
    
    int Ev = i*50;
    stringstream nameCorCor;
    nameCorCor<<"/data/weights/Ev_cor_ME_cor/weights-"<<Ev<<".out";
    stringstream nameCorUnCor;
    nameCorUnCor<<"/data/weights/Ev_cor_ME_uncor/weights-"<<Ev<<".out";
    stringstream nameUnCorCor;
    nameUnCorCor<<"/data/weights/Ev_uncor_ME_cor/weights-"<<Ev<<".out";
    stringstream nameUnCorUnCor;
    nameUnCorUnCor<<"/data/weights/Ev_uncor_ME_uncor/weights-"<<Ev<<".out";
    //Open input file                                                                                                                                          
    ifstream myInfile(nameCorCor.str());
    ifstream myInfile2(nameCorUnCor.str());
    ifstream myInfile3(nameUnCorCor.str());
    ifstream myInfile4(nameUnCorUnCor.str());
    string event_number;
    string likelihood;
    string likelihood_unc;
    float L, L_unc;
    std::vector<float> *Like_corcor = new std::vector<float>();
    std::vector<float> *Like_coruncor = new std::vector<float>();
    std::vector<float> *Like_uncorcor = new std::vector<float>();
    std::vector<float> *Like_uncoruncor = new std::vector<float>();
    
    if ( myInfile.is_open() ) {
      while ( myInfile>>event_number){
        myInfile>>likelihood; myInfile>>likelihood_unc;
        if(atof(likelihood.c_str()) != 0){L = log(atof(likelihood.c_str()));}
        else{L = 999;}
        L_unc = log(atof(likelihood_unc.c_str()));
        Like_corcor->push_back(L);
      }
    }
    myInfile.close();
    
    if ( myInfile2.is_open() ) {
      while ( myInfile2>>event_number){
        myInfile2>>likelihood; myInfile2>>likelihood_unc;
        if(atof(likelihood.c_str()) != 0){L = log(atof(likelihood.c_str()));}
        else{L = 999;}
        L_unc = log(atof(likelihood_unc.c_str()));
        Like_coruncor->push_back(L);
      }
    }
    myInfile2.close();
    
    if( myInfile3.is_open() && myInfile4.is_open() ){
      if ( myInfile3.is_open() ) {
	while ( myInfile3>>event_number){
	  myInfile3>>likelihood; myInfile3>>likelihood_unc;
	  if(atof(likelihood.c_str()) != 0){L = log(atof(likelihood.c_str()));}
	  else{L = 999;}
	  L_unc = log(atof(likelihood_unc.c_str()));
	  Like_uncorcor->push_back(L);
	}
      }
      myInfile3.close();
    
      if ( myInfile4.is_open() ) {
	while ( myInfile4>>event_number){
        myInfile4>>likelihood; myInfile4>>likelihood_unc;
        if(atof(likelihood.c_str()) != 0){L = log(atof(likelihood.c_str()));}
        else{L = 999;}
        L_unc = log(atof(likelihood_unc.c_str()));
        Like_uncoruncor->push_back(L);
	}
      }
      myInfile4.close();
    }
    
  
  for(unsigned int j =0; j < Like_uncorcor->size(); j++){
    if((*Like_corcor)[j] != 999 && (*Like_coruncor)[j] != 999){
      GaussEvCor->Fill((*Like_corcor)[j] - (*Like_coruncor)[j]);
    }
    if((*Like_uncorcor)[j] != 999 && (*Like_uncoruncor)[j] != 999){
      GaussEvUnCor->Fill((*Like_uncorcor)[j] - (*Like_uncoruncor)[j]);
    }
  }
  }
  
  std::vector<float> *Data_uncorcor = new std::vector<float>();
  std::vector<float> *Data_uncoruncor = new std::vector<float>();
  std::vector<float> *Data_cor = new std::vector<float>();
  std::vector<float> *Data_uncor = new std::vector<float>();
  
  
  for(int i=200; i < 8000; i++){
    
    int Ev = i*50;
    stringstream nameCorCor;
    nameCorCor<<"/data/weights/Ev_uncor_ME_cor/weights-"<<Ev<<".out";
    stringstream nameCorUnCor;
    nameCorUnCor<<"/data/weights/Ev_uncor_ME_uncor/weights-"<<Ev<<".out";
    
    //Open input file                                                                                                                                                                                        
    ifstream myInfile(nameCorCor.str());
    ifstream myInfile2(nameCorUnCor.str());
    string event_number;
    string likelihood;
    string likelihood_unc;
    float L, L_unc;
    
    int k = 0;
    if(myInfile.is_open() && myInfile2.is_open()){
      if ( myInfile.is_open() ) {
	while ( myInfile>>event_number){
	  myInfile>>likelihood; myInfile>>likelihood_unc;
	  if(atof(likelihood.c_str()) != 0){L = log(atof(likelihood.c_str()));}
	  else{L = 999;}
	  L_unc = log(atof(likelihood_unc.c_str()));
	  Data_uncorcor->push_back(L);
	  //std::cout<<"event data uncorrelated= "<<Ev+k<<" file number= "<<Ev<<" likelihood correlated= "<<L<<std::endl;                                                                                    
	  k++;
	}
      }
    myInfile.close();
    
    if ( myInfile2.is_open() ) {
      k = 0;
      while ( myInfile2>>event_number){
	myInfile2>>likelihood; myInfile2>>likelihood_unc;
	if(atof(likelihood.c_str()) != 0){L = log(atof(likelihood.c_str()));}
	else{L = 999;}
	L_unc = log(atof(likelihood_unc.c_str()));
	Data_uncoruncor->push_back(L);
	//std::cout<<"event data uncorrelated= "<<Ev+k<<" file number= "<<Ev<<" likelihood uncorrelated= "<<L<<std::endl;                                                                                  
	k++;
      }
    }
    myInfile2.close();
    }
    else{
      for(int x = 0; x < 50; x++)
	{
	  Data_uncorcor->push_back(999);
	  Data_uncoruncor->push_back(999);
	}
    }
  }
  
  for(int i=200; i < 8000; i++){
    
    int Ev = i*50;
    stringstream nameCorCor;
    //   nameCorCor<<"/afs/cern.ch/work/k/ksbeerna/public/weights/Ev_cor_ME_cor/weights-"<<Ev<<".out";
    nameCorCor<<"/data/weights/Ev_cor_ME_cor/weights-"<<Ev<<".out";                                                                                                                                        
    stringstream nameCorUnCor;
    nameCorUnCor<<"/data/weights/Ev_cor_ME_uncor/weights-"<<Ev<<".out";
    
    //Open input file                                                                                                                                                                                        
    ifstream myInfile(nameCorCor.str());
    ifstream myInfile2(nameCorUnCor.str());
    string event_number;
    string likelihood;
    string likelihood_unc;
    float L, L_unc;
    int k = 0;
    if(myInfile.is_open() && myInfile2.is_open()){
      if ( myInfile.is_open() ) {
	while ( myInfile>>event_number){
	  myInfile>>likelihood; myInfile>>likelihood_unc;
	  if(atof(likelihood.c_str()) != 0){L = log(atof(likelihood.c_str()));}
	  else{L = 999;}
	  L_unc = log(atof(likelihood_unc.c_str()));
	  Data_cor->push_back(L);
	  //std::cout<<"event data correlated= "<<Ev+k<<" file number= "<<Ev<<" likelihood correlated= "<<L<<std::endl;                                                                                      
	  k++;
	}
      }
      myInfile.close();
      
      k = 0;
      if ( myInfile2.is_open() ) {
	while ( myInfile2>>event_number){
	  myInfile2>>likelihood; myInfile2>>likelihood_unc;
	  if(atof(likelihood.c_str()) != 0){L = log(atof(likelihood.c_str()));}
	  else{L = 999;}
	  L_unc = log(atof(likelihood_unc.c_str()));
	  Data_uncor->push_back(L);
	  //std::cout<<"event data correlated= "<<Ev+k<<" file number= "<<Ev<<" likelihood uncorrelated= "<<L<<std::endl;                                                                                      
	  k++;
	}
      }
      myInfile2.close();
    }
    else{
      for(int x = 0; x < 50; x++)
	{
	  Data_cor->push_back(999);
	  Data_uncor->push_back(999);
	}
    }
  }
  
  TH1F* Pull = new TH1F("Pull","(f_obs - f_exp)/sigma_obs", 50, -3, 3);
   std::vector<float> *Lev = new std::vector<float>();
  std::vector<TH1F> *hists = new std::vector<TH1F>();
  TH1F* Data = new TH1F("Data","-2ln#lambda = 2lnL_{H=C} - 2lnL_{H=U}", 100, -3, 3);
  hists->push_back(*Data);
  int i = 0;
  int x = 0;
 
  for(unsigned int j =0; j < Data_cor->size(); j++){
    if((*Data_cor)[j] != 999 && (*Data_uncor)[j] != 999){Lev->push_back((*Data_cor)[j] - (*Data_uncor)[j]);                   
    (*hists)[i].Fill((*Lev)[x]);
    x++;}
    if((x + 1) % 2000 == 0){
      if((*hists)[i].GetRMS() >= 0.5 ){std::cout<<i*50<<std::endl;}
      i++;
      TH1F* Data2 = new TH1F("Data2","-2ln#lambda = 2lnL_{H=C} - 2lnL_{H=U}", 100, -3, 3);
      hists->push_back(*Data2);
    }
  }
  
  for(unsigned int z = 0; z < hists->size(); z++){
    TCanvas *c1 = new TCanvas();
    TH1F H = (*hists)[z];
    H.Sumw2();
    H.Scale(1./H.Integral("width"));
    std::vector<double> res = DoFit(&H, GaussEvCor, GaussEvUnCor);
    double f= (res)[0];
    double Err_f=(res)[1];
    Pull->Fill((f - 1.0)/Err_f);
    std::cout<<"f= "<<f<<" error= "<<Err_f<<std::endl;
    H.Draw();
    c1->Modified();
    c1->Update();
    stringstream namez;
    namez<<"PullTemplate"<<z<<".png";
    c1->SaveAs((TString) namez.str());
    delete c1;
  }
  
  TCanvas *c2 = new TCanvas();
  Pull->Draw();
  c2->Modified();
  c2->Update();
  c2->SaveAs("PullDataTemplateCorrelated.C");
  delete c2;

  return 1;
}


