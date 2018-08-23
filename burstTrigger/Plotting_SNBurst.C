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
std::vector<int> SelectedConfig = {0,1,2,3,4,5,10};

int main()
{
  TString s_Filename = "GH_SNMC";
  ifstream inFile;
  inFile.open("Analyse_SNBurst_"+s_Filename+".txt");

  std::map<std::pair<int, double>,std::pair<double,double>> map_ConfigToEffAndBkgd;
  int Config;
  double TimeWindow, Eff, Bkgd, dummy;
  std::vector<int> vec_Config;

  while(inFile >> Config >> TimeWindow>> Eff >> Bkgd >> dummy){
  std::cout << "Config " << Config << " TimeWindow " << TimeWindow
            << ", Eff " << Eff
            << ", Bkgd " << Bkgd
            << std::endl;

    if(std::find(vec_Config.begin(), vec_Config.end(), Config) == vec_Config.end())
      vec_Config.push_back(Config);

    map_ConfigToEffAndBkgd[std::make_pair(Config,TimeWindow)] = {Eff,Bkgd};
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
  std::map<std::pair<int,double>,TH1D*>   map_h_FakeRateVNClusters;
  std::map<std::pair<int,double>,TH1D*>   map_h_FakeRateVNClustersLow;
  std::map<std::pair<int,double>,TH1D*>   map_h_EfficiencyVEvents;
  std::map<std::pair<int,double>,TH1D*>   map_h_EfficiencyVDistance;
  std::map<std::pair<int,double>,TH1D*>   map_h_EffGalaxy;
  std::map<std::pair<int,double>,TGraph*> map_g_ROC;

  std::map<int,std::pair<double,double>>::iterator it_Config;
  int globalIt=0;
  for(auto& it_Config : map_ConfigToEffAndBkgd){

    int color = vec_Colors.at(globalIt % vec_Colors.size());

    int config     = it_Config.first.first;
    TimeWindow     = it_Config.first.second;
    std::cout << "Dealing with Config " << config <<  " time window " << TimeWindow << std::endl;
    TString s_FakeRateVNClusters     = Form("h_FakeRateVNClusters_Config%i_TW%.3f",    config, TimeWindow);
    TString s_FakeRateVNClustersLow  = Form("h_FakeRateVNClustersLow_Config%i_TW%.3f", config, TimeWindow);
    TString s_EfficiencyVEvents      = Form("h_EfficiencyVEvents_Config%i_TW%.3f",     config, TimeWindow);
    TString s_EfficiencyVDistance    = Form("h_EfficiencyVDistance_Config%i_TW%.3f",   config, TimeWindow);
    TString s_EffGalaxy              = Form("h_NeighbourhoodEffiency_Config%i_TW%.3f", config, TimeWindow);
    TString s_ROC                    = Form("g_ROC_Config%i_TW%.3f",                   config, TimeWindow);

    TH1D *h_FakeRateVNClusters    = (TH1D*)f_Input->Get(s_FakeRateVNClusters);
    TH1D *h_FakeRateVNClustersLow = (TH1D*)f_Input->Get(s_FakeRateVNClustersLow);

    TH1D *h_EfficiencyVEvents = (TH1D*)f_Input->Get(s_EfficiencyVEvents);
    if(!h_EfficiencyVEvents){
      std::cout << "Erasing config " << config << " as there is no EfficiencyVEvents plot" << std::endl;
      map_ConfigToEffAndBkgd.erase(std::make_pair(config,TimeWindow));
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
    if(!h_EffGalaxy){
      std::cout << "Erasing config " << config << " as there is no NeighbourhoodEffiency plot" << std::endl;
      map_ConfigToEffAndBkgd.erase(std::make_pair(config,TimeWindow));
      continue;
    }
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

    if (config>=100) {
      h_FakeRateVNClusters   ->SetLineStyle(2);
      h_FakeRateVNClustersLow->SetLineStyle(2);
      h_EfficiencyVEvents    ->SetLineStyle(2);
      h_EfficiencyVDistance  ->SetLineStyle(2);
      h_EffGalaxy            ->SetLineStyle(2);
      g_ROC                  ->SetLineStyle(2);
    }
    g_ROC->SetMarkerStyle(3);

    map_h_FakeRateVNClusters   [std::make_pair(config,TimeWindow)] = h_FakeRateVNClusters;
    map_h_FakeRateVNClustersLow[std::make_pair(config,TimeWindow)] = h_FakeRateVNClustersLow;
    map_h_EfficiencyVEvents    [std::make_pair(config,TimeWindow)] = h_EfficiencyVEvents;
    map_h_EfficiencyVDistance  [std::make_pair(config,TimeWindow)] = h_EfficiencyVDistance;
    map_h_EffGalaxy            [std::make_pair(config,TimeWindow)] = h_EffGalaxy;
    map_g_ROC                  [std::make_pair(config,TimeWindow)] = g_ROC;
    globalIt++;
  }
  TCanvas *c_Global = new TCanvas();
  c_Global->Print("Results.pdf[");
  std::string legHeader = "Individual Marley Eff & 10kt Bkgd Rate";
  std::string legEntryFormatWire = "Wire Clusters - Eff: %.2f & Bkgd rate: %.2f Hz (%i s timing window)";
  std::string legEntryFormatOpti = "Optical Custers (nHit>= %i): - Eff: %.2f & Bkgd rate: %.2f Hz (%i s timing window)";

  THStack *stk_FakeRateVNClusters = new THStack("stk_FakeRateVNClusters", "Number of Clusters in Time Window Required to Trigger vs. Trigger Rate");
  TLegend *leg_FakeRateVNClusters = new TLegend(0.05, 0.05, 0.95, 0.95);
  leg_FakeRateVNClusters->SetHeader(legHeader.c_str());

  for(auto const& it_Config : map_ConfigToEffAndBkgd){
    stk_FakeRateVNClusters->Add(map_h_FakeRateVNClusters     [it_Config.first]);
    stk_FakeRateVNClusters->Add(map_h_FakeRateVNClustersLow  [it_Config.first]);
    if(it_Config.first.first>=100) {
      leg_FakeRateVNClusters->AddEntry(map_h_FakeRateVNClusters[it_Config.first],
                                       Form(legEntryFormatWire.c_str(), it_Config.second.first, it_Config.second.second, (int)it_Config.first.second), "L");
    } else {
      leg_FakeRateVNClusters->AddEntry(map_h_FakeRateVNClusters[it_Config.first],
                                       Form(legEntryFormatOpti.c_str(), it_Config.first.first, it_Config.second.first, it_Config.second.second, (int)it_Config.first.second), "L");
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
  TText *t_perYear = new TText(map_h_FakeRateVNClusters.begin()->second->GetXaxis()->GetXmax()-0.15*range, 1.7e-8,    "1/Year");
  TText *t_perMonth  = new TText(map_h_FakeRateVNClusters.begin()->second->GetXaxis()->GetXmax()-0.15*range, 2e-7,    "1/Month");
  TLine *l_perYear = new TLine(0, 3.44e-8, map_h_FakeRateVNClusters.begin()->second->GetXaxis()->GetXmax(), 3.44e-8);
  TLine *l_perMonth  = new TLine(0, 4.13e-7, map_h_FakeRateVNClusters.begin()->second->GetXaxis()->GetXmax(), 4.13e-7);

  l_perYear->SetLineColor(4);
  l_perMonth->SetLineColor(4);
  l_perYear->SetLineWidth(3);
  l_perMonth->SetLineWidth(3);
  t_perYear->Draw();
  t_perMonth->Draw();
  l_perYear->Draw();
  l_perMonth->Draw();
  gPad->RedrawAxis();
  c_Global->Print("Results.pdf");

  c_Global->SetLogx(false);
  THStack *stk_EfficiencyVEvents = new THStack("stk_EfficiencyVEvents", "Efficiency vs. Number of Events in SN Burst, Fake Trigger Rate: 1/Year");
  for(auto const& it_Config : map_ConfigToEffAndBkgd){
    stk_EfficiencyVEvents->Add(map_h_EfficiencyVEvents[it_Config.first]);
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

  THStack *stk_EfficiencyVDistance = new THStack("stk_EfficiencyVDistance", "Efficiency vs. Distance to SN, Fake Trigger Rate: 1/Year");
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
  THStack *stk_EffGalaxy = new THStack("stk_EffGalaxy", "Galactic Neighbourhood Coverage, Fake Trigger Rate 1/Year");

  stk_EffGalaxy->Add(h_SNProbabilityVDistance);

  for(auto const& it_Config : map_ConfigToEffAndBkgd){
    stk_EffGalaxy->Add(map_h_EffGalaxy[it_Config.first]);
  }
  c_Global->Clear();
  c_Global->Draw();
  c_Global->SetLogx(false);
  c_Global->SetLogy();
  stk_EffGalaxy->SetMinimum(1e-5);
  stk_EffGalaxy->Draw("NOSTACK");
  stk_EffGalaxy->GetXaxis()->SetTitle("SN Distance, (kpc)");
  stk_EffGalaxy->GetYaxis()->SetTitle("Efficiency x SN Probability");
  stk_EffGalaxy->Draw("NOSTACK HIST");
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
  TLine *l_perYear_2 = new TLine(minX, 3.44e-8, maxX, 3.44e-8);
  l_perYear_2->SetLineColor(1);
  l_perYear_2->SetLineWidth(3);
  TText *t_perYear_2 = new TText(minX+0.015, 7e-8, "1/Year");
  gPad->RedrawAxis();
  l_perYear_2->Draw();
  t_perYear_2->Draw();
  c_Global->Print("Results.pdf");
  c_Global->Print("Results.pdf]");

  return 0;
}
