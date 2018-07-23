#include <algorithm>
#include <iostream>
#include <vector>
#include <fstream>

#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TText.h"
#include "TROOT.h"

#include "colors.h"

using namespace std;
std::vector<int> SelectedConfig = {0,1,2,3,4,5};

int main()
{
  TString s_Filename = "GH_SNMC";
  ifstream inFile;
  inFile.open("Analyse_SNBurst_"+s_Filename+".txt");

  std::map<int,std::pair<double,double>> map_ConfigToEffAndBkgd;
  int Config;
  double Eff, Bkgd, dummy;
  std::vector<int> vec_Config;
  
  while(inFile >> Config >> Eff >> Bkgd >> dummy){
    std::cout << "Config " << Config
      << ", Eff " << Eff
      << ", Bkgd " << Bkgd
      << std::endl;

    if(std::find(vec_Config.begin(), vec_Config.end(), Config) == vec_Config.end())
      vec_Config.push_back(Config);
    
    map_ConfigToEffAndBkgd[Config] = {Eff,Bkgd};
  }
  
  int nConfig = (int)vec_Config.size();
  std::cout << "There are " << nConfig << " configs." << std::endl;
  
  gROOT->ForceStyle();
  std::vector<int> vec_Colors = getColors(2);
  TFile *f_SNTheoryDistributions = new TFile("SNTheoryDistributions.root","READ");
  TH1D  *h_SNProbabilityVDistance = (TH1D*)f_SNTheoryDistributions->Get("h_SNProbabilityVDistance_LMC");
  h_SNProbabilityVDistance->SetLineWidth(3);
  h_SNProbabilityVDistance->SetLineColor(46);

  TFile *f_Input  = new TFile("Analyse_SNBurst_GH_SNMC.root", "READ");
  std::map<int,TH1D*>   map_h_FakeRateVNClusters;
  std::map<int,TH1D*>   map_h_FakeRateVNClustersLow;
  std::map<int,TH1D*>   map_h_EfficiencyVEvents;
  std::map<int,TH1D*>   map_h_EfficiencyVDistance;
  std::map<int,TH1D*>   map_h_EffGalaxy;
  std::map<int,TGraph*> map_g_ROC;

  std::map<int,std::pair<double,double>>::iterator it_Config;
  int globalIt=0;
  for(auto const& it_Config : map_ConfigToEffAndBkgd){
    int color = vec_Colors.at(globalIt % vec_Colors.size());

    int config     = it_Config.first;
    
    std::cout << "Dealing with Config " << config << std::endl;
    TString s_FakeRateVNClusters     = Form("h_FakeRateVNClusters_Config%i",    config);
    TString s_FakeRateVNClustersLow  = Form("h_FakeRateVNClustersLow_Config%i", config);
    TString s_EfficiencyVEvents      = Form("h_EfficiencyVEvents_Config%i",     config);
    TString s_EfficiencyVDistance    = Form("h_EfficiencyVDistance_Config%i",   config);
    TString s_EffGalaxy              = Form("h_NeighbourhoodEffiency_Config%i", config);
    TString s_ROC                    = Form("g_ROC_Config%i",                   config);
    
    TH1D *h_FakeRateVNClusters    = (TH1D*)f_Input->Get(s_FakeRateVNClusters); 
    TH1D *h_FakeRateVNClustersLow = (TH1D*)f_Input->Get(s_FakeRateVNClustersLow);
    
    TH1D *h_EfficiencyVEvents = (TH1D*)f_Input->Get(s_EfficiencyVEvents);
    if(!h_EfficiencyVEvents){
      std::cout << "Erasing config " << config << " as there is no EfficiencyVEvents plot" << std::endl;
      map_ConfigToEffAndBkgd.erase(config);
      continue;
    }

    int   maxBin_Events = h_EfficiencyVEvents->GetMaximumBin();
    for(int i = maxBin_Events; i < h_EfficiencyVEvents->GetSize()-1; i++)
      h_EfficiencyVEvents->SetBinContent(i,1);

    TH1D *h_EfficiencyVDistance = (TH1D*)f_Input->Get(s_EfficiencyVDistance);
    int   maxBin_Distance = h_EfficiencyVDistance->GetMaximumBin();
    for(int i = maxBin_Distance; i > 0; i--)
      h_EfficiencyVDistance->SetBinContent(i,1);
    TH1D   *h_EffGalaxy           = (TH1D*)f_Input->Get(s_EffGalaxy);
    TGraph *g_ROC                 = (TGraph*)f_Input->Get(s_ROC);
    h_FakeRateVNClusters   ->SetLineColor  (color);
    h_FakeRateVNClustersLow->SetLineColor  (color);
    h_FakeRateVNClusters   ->SetMarkerColor(color);
    h_FakeRateVNClustersLow->SetMarkerColor(color);
    h_EfficiencyVEvents    ->SetLineColor  (color);
    h_EfficiencyVDistance  ->SetLineColor  (color);
    h_EffGalaxy            ->SetLineColor  (color);
    g_ROC                  ->SetMarkerColor(color);
    g_ROC                  ->SetLineColor  (color);

    g_ROC->SetMarkerStyle(3);

    map_h_FakeRateVNClusters   [config] = h_FakeRateVNClusters;
    map_h_FakeRateVNClustersLow[config] = h_FakeRateVNClustersLow;
    map_h_EfficiencyVEvents    [config] = h_EfficiencyVEvents;
    map_h_EfficiencyVDistance  [config] = h_EfficiencyVDistance;
    map_h_EffGalaxy            [config] = h_EffGalaxy;
    map_g_ROC                  [config] = g_ROC;
    globalIt++;    
  }
  TCanvas *c_Global = new TCanvas();
  c_Global->Print("Results.pdf[");
  std::string legHeader = "Individual Marley Eff & 10kt Bkgd Rate";
  std::string legEntryFormatWire = "Wire Clusters - Eff: %.2f & Bkgd rate: %.2f Hz 10s timing windown";
  std::string legEntryFormatOpti = "Optical Custers (nHit>= %i): - Eff: %.2f & Bkgd rate: %.2f Hz 1s timing windown";
  
  THStack *stk_FakeRateVNClusters = new THStack("stk_FakeRateVNClusters", "Number of Clusters in Time Window Required to Trigger vs. Trigger Rate");
  TLegend *leg_FakeRateVNClusters = new TLegend(0.05, 0.05, 0.95, 0.95);
  leg_FakeRateVNClusters->SetHeader(legHeader.c_str());

  for(auto const& it_Config : map_ConfigToEffAndBkgd){
    stk_FakeRateVNClusters->Add(map_h_FakeRateVNClusters     [it_Config.first]);
    stk_FakeRateVNClusters->Add(map_h_FakeRateVNClustersLow  [it_Config.first]);
    if(it_Config.first==0) {
      leg_FakeRateVNClusters->AddEntry(map_h_FakeRateVNClusters[it_Config.first],
                                       Form(legEntryFormatWire.c_str(), it_Config.second.first, it_Config.second.second), "L");
    } else {
      leg_FakeRateVNClusters->AddEntry(map_h_FakeRateVNClusters[it_Config.first],
                                       Form(legEntryFormatOpti.c_str(), it_Config.first, it_Config.second.first, it_Config.second.second), "L");
    }
  }
  leg_FakeRateVNClusters->Draw();
  c_Global->Print("Results.pdf");

  c_Global->SetLogy();
  c_Global->SetLogx();
  stk_FakeRateVNClusters->SetMinimum(1e-9);
  stk_FakeRateVNClusters->Draw("NOSTACK C");
  
  stk_FakeRateVNClusters->GetXaxis()->SetTitle("Number of Clusters/Time Window");
  stk_FakeRateVNClusters->GetYaxis()->SetTitle("Trigger Rate, (Hz)");
  stk_FakeRateVNClusters->GetXaxis()->SetLimits(0,map_h_FakeRateVNClusters.begin()->second->GetXaxis()->GetXmax());
  double range = map_h_FakeRateVNClusters.begin()->second->GetXaxis()->GetXmax() - map_h_FakeRateVNClusters.begin()->second->GetXaxis()->GetXmin();
  TText *t_perMonth  = new TText(map_h_FakeRateVNClusters.begin()->second->GetXaxis()->GetXmax()-0.15*range, 2e-7,    "1/Month");
  TText *t_perWeek   = new TText(map_h_FakeRateVNClusters.begin()->second->GetXaxis()->GetXmax()-0.15*range, 1.65e-6, "1/Week");
  TText *t_perDay    = new TText(map_h_FakeRateVNClusters.begin()->second->GetXaxis()->GetXmax()-0.15*range, 1.16e-5, "1/Day");
  TLine *l_perMonth  = new TLine(0, 4.13e-7, map_h_FakeRateVNClusters.begin()->second->GetXaxis()->GetXmax(), 4.13e-7);
  TLine *l_perWeek   = new TLine(0, 1.65e-6, map_h_FakeRateVNClusters.begin()->second->GetXaxis()->GetXmax(), 1.65e-6);
  TLine *l_perDay    = new TLine(0, 1.16e-5, map_h_FakeRateVNClusters.begin()->second->GetXaxis()->GetXmax(), 1.16e-5);
  l_perMonth->SetLineColor(4);
  l_perWeek ->SetLineColor(4);
  l_perDay  ->SetLineColor(4);
  l_perMonth->SetLineWidth(3);
  l_perWeek ->SetLineWidth(3);
  l_perDay  ->SetLineWidth(3);
  t_perMonth->Draw();
  t_perWeek ->Draw();
  t_perDay  ->Draw();
  l_perMonth->Draw();
  l_perWeek ->Draw();
  l_perDay  ->Draw();
  gPad->RedrawAxis(); 
  c_Global->Print("Results.pdf");

  c_Global->SetLogx(false);
  THStack *stk_EfficiencyVEvents = new THStack("stk_EfficiencyVEvents", "Efficiency vs. Number of Events in SN Burst, Fake Trigger Rate: 1/Month");
  for(auto const& it_Config : map_ConfigToEffAndBkgd){
    stk_EfficiencyVEvents->Add(map_h_EfficiencyVEvents     [it_Config.first]);
  }
  c_Global->Clear();
  c_Global->SetLogy(false);
  c_Global->Draw();
  c_Global->SetLogx();
  stk_EfficiencyVEvents->Draw("NOSTACK");
  stk_EfficiencyVEvents->GetXaxis()->SetTitle("Number of Events in SN Burst");
  stk_EfficiencyVEvents->GetYaxis()->SetTitle("Efficiency");
  gPad->RedrawAxis();
  c_Global->Print("Results.pdf");

  THStack *stk_EfficiencyVDistance = new THStack("stk_EfficiencyVDistance", "Efficiency vs. Distance to SN, Fake Trigger Rate: 1/Month");
  for(auto const& it_Config : map_ConfigToEffAndBkgd){
    stk_EfficiencyVDistance->Add(map_h_EfficiencyVDistance[it_Config.first]);
  }
  c_Global->Clear();
  c_Global->Draw();
  c_Global->SetLogx();
  stk_EfficiencyVDistance->Draw("NOSTACK");
  stk_EfficiencyVDistance->GetXaxis()->SetTitle("SN Distance, (kpc)");
  stk_EfficiencyVDistance->GetYaxis()->SetTitle("Efficiency");
  gPad->RedrawAxis();
  c_Global->Print("Results.pdf");


  gStyle->SetOptStat(0);
  THStack *stk_EffGalaxy = new THStack("stk_EffGalaxy", "Galactic Neighbourhood Coverage, Fake Trigger Rate 1/Month");

  stk_EffGalaxy->Add(h_SNProbabilityVDistance);

  for(auto const& it_Config : map_ConfigToEffAndBkgd){
    stk_EffGalaxy->Add(map_h_EffGalaxy[it_Config.first]);
  }
  c_Global->Clear();
  c_Global->Draw();
  c_Global->SetLogx(false);
  c_Global->SetLogy();
  stk_EffGalaxy->Draw("NOSTACK");
  stk_EffGalaxy->GetXaxis()->SetTitle("SN Distance, (kpc)");
  stk_EffGalaxy->GetYaxis()->SetTitle("Efficiency x SN Probability");
  stk_EffGalaxy->Draw("NOSTACK");
  gPad->RedrawAxis();
  c_Global->Print("Results.pdf");


  c_Global->Clear();
  c_Global->Draw();
  double minX=0.6, maxX=1;
  double minY=10e-15, maxY=10;

  map_g_ROC.begin()->second->GetXaxis()->SetLimits(minX, maxX);
  map_g_ROC.begin()->second->SetMaximum(maxY);
  map_g_ROC.begin()->second->SetMinimum(minY);
  map_g_ROC.begin()->second->SetTitle("Fake Trigger Rate vs. Galactic Neighbourhood Coverage");
  map_g_ROC.begin()->second->GetXaxis()->SetTitle("Galactic Neighbourhood Coverage");
  map_g_ROC.begin()->second->GetYaxis()->SetTitle("Fake Trigger Rate, (Hz)");
  map_g_ROC.begin()->second->SetMaximum(10);
  map_g_ROC.begin()->second->Draw("AP");

  for(auto const& it_Config : map_ConfigToEffAndBkgd){
    map_g_ROC[it_Config.first]->Draw("P");
  }
  TLine *l_perMonth_2 = new TLine(minX, 4.13e-7, maxX, 4.13e-7);
  l_perMonth_2->SetLineColor(1);
  l_perMonth_2->SetLineWidth(3);
  TText *t_perMonth_2 = new TText(minX+0.015, 8e-7, "1/Month");
  gPad->RedrawAxis();
  l_perMonth_2->Draw();
  t_perMonth_2->Draw();
  c_Global->Print("Results.pdf");
  c_Global->Print("Results.pdf]");

  return 0;
}
