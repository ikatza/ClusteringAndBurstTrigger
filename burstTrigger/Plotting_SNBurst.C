#include <algorithm>
#include <iostream>
#include <vector>
#include <fstream>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TStyle.h>
#include <TText.h>
#include <TROOT.h>
#include "/nashome/p/plasorak/utils/colors.h"

using namespace std;
std::vector<int> SelectedConfig = {3};


int main()
{
  TString s_Filename = "GH_SNMC";
  ifstream inFile;
  inFile.open("Analyse_SNBurst_"+s_Filename+".txt");

  TString s_Filename2 = "ALEX";
  ifstream inFile2;
  inFile2.open("Analyse_SNBurst_"+s_Filename2+".txt");

  std::map<std::pair<int,int>,std::pair<double,double>> map_ConfigToEffAndBkgd;
  int TimeWindow;
  int Config;
  double Eff, Bkgd, dummy;
  std::vector<double> vec_TimeWindow;
  std::vector<int> vec_Config;
  
  while(inFile >> TimeWindow >> Config >> Eff >> Bkgd >> dummy){
    std::cout << "file1m TimeWindow " << TimeWindow
      << ", Config " << Config
      << ", Eff " << Eff
      << ", Bkgd " << Bkgd
      << std::endl;
    if(TimeWindow <= 0) continue; // ONLY TREAT THE GOOD TW

    // if(std::find(SelectedConfig.begin(), SelectedConfig.end(), Config) == SelectedConfig.end())
    //   continue; // ONLY TREAT SELECTED CONFIGURATIONS

    if(std::find(vec_TimeWindow.begin(), vec_TimeWindow.end(), TimeWindow) == vec_TimeWindow.end())
      vec_TimeWindow.push_back(TimeWindow);

    if(std::find(vec_Config.begin(), vec_Config.end(), Config) == vec_Config.end())
      vec_Config.push_back(Config);
    
    map_ConfigToEffAndBkgd[{TimeWindow,Config}] = {Eff,Bkgd};
  }
  
  std::map<std::pair<int,int>,std::pair<double,double>> map_ConfigToEffAndBkgd2;
  while(inFile2 >> TimeWindow >> Config >> Eff >> Bkgd >> dummy){
    std::cout << "file2 TimeWindow " << TimeWindow
      << ", Config " << Config
      << ", Eff " << Eff
      << ", Bkgd " << Bkgd
      << std::endl;
    if(TimeWindow <= 0) continue; // ONLY TREAT THE GOOD TW
    map_ConfigToEffAndBkgd2[{TimeWindow,Config}] = {Eff,Bkgd};
  }

  int nTimeWindow = (int)vec_TimeWindow.size();
  std::cout << "There are " << nTimeWindow << " time windows." << std::endl;
  
  int nConfig = (int)vec_Config.size();
  std::cout << "There are " << nConfig << " configs." << std::endl;
  
  gROOT->ForceStyle();
  std::vector<int> vec_Colors = getColors(2);
  TFile *f_SNTheoryDistributions = new TFile("SNTheoryDistributions.root","READ");
  TH1D  *h_SNProbabilityVDistance = (TH1D*)f_SNTheoryDistributions->Get("h_SNProbabilityVDistance_LMC");
  h_SNProbabilityVDistance->SetLineWidth(3);
  h_SNProbabilityVDistance->SetLineColor(46);

  TFile *f_Input  = new TFile("Analyse_SNBurst_GH_SNMC.root", "READ");
  TFile *f_Input2 = new TFile("Analyse_SNBurst_ALEX.root", "READ");
  std::cout << "Alex's config" << std::endl;
  std::cout << "config " << map_ConfigToEffAndBkgd2.begin()->first.second << std::endl;
  std::cout << "tw     " << map_ConfigToEffAndBkgd2.begin()->first.first << std::endl;
  
  TH1D*   h_Alex_FakeRateVNClusters    = (TH1D*  )f_Input2->Get(Form("h_FakeRateVNClusters_Config%i_TimeWindow%i",
                                                                     map_ConfigToEffAndBkgd2.begin()->first.second,
                                                                     map_ConfigToEffAndBkgd2.begin()->first.first));
  TH1D*   h_Alex_FakeRateVNClustersLow = (TH1D*  )f_Input2->Get(Form("h_FakeRateVNClustersLow_Config%i_TimeWindow%i",
                                                                     map_ConfigToEffAndBkgd2.begin()->first.second,
                                                                     map_ConfigToEffAndBkgd2.begin()->first.first));
  TH1D*   h_Alex_EfficiencyVEvents     = (TH1D*  )f_Input2->Get(Form("h_EfficiencyVEvents_Config%i_TimeWindow%i",
                                                                     map_ConfigToEffAndBkgd2.begin()->first.second,
                                                                     map_ConfigToEffAndBkgd2.begin()->first.first));
  TH1D*   h_Alex_EfficiencyVDistance   = (TH1D*  )f_Input2->Get(Form("h_EfficiencyVDistance_Config%i_TimeWindow%i",
                                                                     map_ConfigToEffAndBkgd2.begin()->first.second,
                                                                     map_ConfigToEffAndBkgd2.begin()->first.first));
  TH1D*   h_Alex_EffGalaxy             = (TH1D*  )f_Input2->Get(Form("h_NeighbourhoodEffiency_Config%i_TimeWindow%i",
                                                                     map_ConfigToEffAndBkgd2.begin()->first.second,
                                                                     map_ConfigToEffAndBkgd2.begin()->first.first));
  TGraph* g_Alex_ROC                   = (TGraph*)f_Input2->Get(Form("g_ROC_Config%i_TimeWindow%i",
                                                                     map_ConfigToEffAndBkgd2.begin()->first.second,
                                                                     map_ConfigToEffAndBkgd2.begin()->first.first));

  

  
  std::map<std::pair<int,int>,TH1D*>   map_h_FakeRateVNClusters;
  std::map<std::pair<int,int>,TH1D*>   map_h_FakeRateVNClustersLow;
  std::map<std::pair<int,int>,TH1D*>   map_h_EfficiencyVEvents;
  std::map<std::pair<int,int>,TH1D*>   map_h_EfficiencyVDistance;
  std::map<std::pair<int,int>,TH1D*>   map_h_EffGalaxy;
  std::map<std::pair<int,int>,TGraph*> map_g_ROC;
  std::map<int, TGraph*>               map_g_EfficiencyVTimeWindow;

  std::vector<int>::iterator it_Config2;

  for(it_Config2  = vec_Config.begin();
      it_Config2 != vec_Config.end(); ++it_Config2){
    int config=*it_Config2;
    std::cout << "Creating " << config << std::endl;
    map_g_EfficiencyVTimeWindow[config] = new TGraph(vec_TimeWindow.size()-2);
    map_g_EfficiencyVTimeWindow[config]->SetName(Form("g_EffVTime_Config%i", config));
    map_g_EfficiencyVTimeWindow[config]->SetMarkerColor(vec_Colors[config]);
    map_g_EfficiencyVTimeWindow[config]->SetMarkerStyle(3);
    map_g_EfficiencyVTimeWindow[config]->SetLineColor  (vec_Colors[config]);
  }

  std::map<std::pair<int,int>,std::pair<double,double>>::iterator it_Config;
  int globalIt=0;
  for(it_Config  = map_ConfigToEffAndBkgd.begin();
      it_Config != map_ConfigToEffAndBkgd.end(); ++it_Config){
    int color = vec_Colors.at(globalIt % vec_Colors.size());

    int config     = (*it_Config).first.second;
    int TimeWindow = (*it_Config).first.first;
  
    std::cout << "Dealing with Config " << config << ", time window: " << TimeWindow << std::endl;
    TString s_FakeRateVNClusters     = Form("h_FakeRateVNClusters_Config%i_TimeWindow%i",    config, TimeWindow);
    TString s_FakeRateVNClustersLow  = Form("h_FakeRateVNClustersLow_Config%i_TimeWindow%i", config, TimeWindow);
    TString s_EfficiencyVEvents      = Form("h_EfficiencyVEvents_Config%i_TimeWindow%i",     config, TimeWindow);
    TString s_EfficiencyVDistance    = Form("h_EfficiencyVDistance_Config%i_TimeWindow%i",   config, TimeWindow);
    TString s_EffGalaxy              = Form("h_NeighbourhoodEffiency_Config%i_TimeWindow%i", config, TimeWindow);
    TString s_ROC                    = Form("g_ROC_Config%i_TimeWindow%i",                   config, TimeWindow);
    
    TH1D *h_FakeRateVNClusters    = (TH1D*)f_Input->Get(s_FakeRateVNClusters); 
    TH1D *h_FakeRateVNClustersLow = (TH1D*)f_Input->Get(s_FakeRateVNClustersLow); 
    
    TH1D *h_EfficiencyVEvents = (TH1D*)f_Input->Get(s_EfficiencyVEvents);
    if(!h_EfficiencyVEvents){
      map_ConfigToEffAndBkgd.erase(it_Config);
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

    map_h_FakeRateVNClusters   [{TimeWindow,config}] = h_FakeRateVNClusters;
    map_h_FakeRateVNClustersLow[{TimeWindow,config}] = h_FakeRateVNClustersLow;
    map_h_EfficiencyVEvents    [{TimeWindow,config}] = h_EfficiencyVEvents;
    map_h_EfficiencyVDistance  [{TimeWindow,config}] = h_EfficiencyVDistance;
    map_h_EffGalaxy            [{TimeWindow,config}] = h_EffGalaxy;
    map_g_ROC                  [{TimeWindow,config}] = g_ROC;
    globalIt++;    
  }
  TCanvas *c_Global = new TCanvas();
  c_Global->Print("Results.pdf[");
  std::vector<double>::iterator it_TimeWindow2;
  it_Config2;
  int j=0;
  for(it_Config2  = vec_Config.begin();
      it_Config2 != vec_Config.end(); ++it_Config2){
    map_g_EfficiencyVTimeWindow[*it_Config2]->SetMarkerColor(vec_Colors[j]);
    j++;
    int i=0;
    for(it_TimeWindow2  = vec_TimeWindow.begin();
        it_TimeWindow2 != vec_TimeWindow.end(); ++it_TimeWindow2){
      if((*it_TimeWindow2)<=100)continue;
      TH1D* eff = map_h_EfficiencyVDistance[{*it_TimeWindow2,*it_Config2}];
      int bin = eff->FindBin(50.);
      double binc = eff->GetBinContent(bin);
      std::cout << " binc " << binc << std::endl;
      map_g_EfficiencyVTimeWindow[*it_Config2]->SetPoint(i, *it_TimeWindow2,binc);
      i++;
    }
  }

  std::string legHeader = "Individual Marley Eff & 10kt Bkgd Rate";
  std::string legEntryFormat = "Eff: %.2f & Bkgd rate: %.2f Hz";
  
  map_g_EfficiencyVTimeWindow[0]->SetMaximum(5);
  gPad->SetLogy();
  //gPad->SetLogx();
  map_g_EfficiencyVTimeWindow[0]->SetMinimum(0.01);
  map_g_EfficiencyVTimeWindow[0]->Draw("AP");
  map_g_EfficiencyVTimeWindow[0]->SetTitle("");
  map_g_EfficiencyVTimeWindow[0]->GetXaxis()->SetTitle("Time window [ms]");
  map_g_EfficiencyVTimeWindow[0]->GetYaxis()->SetTitle("Efficiency");
  TLegend *leg_EfficiencyVTimeWindow = new TLegend(0.5, 0.7, 0.9, 0.9);
  for(it_Config2  = vec_Config.begin();
      it_Config2 != vec_Config.end(); ++it_Config2){
    map_g_EfficiencyVTimeWindow[*it_Config2]->Draw("PC");
    std::cout << map_ConfigToEffAndBkgd[{5000,*it_Config2}].second << std::endl;
    leg_EfficiencyVTimeWindow->AddEntry(map_g_EfficiencyVTimeWindow[*it_Config2],
                                        Form(legEntryFormat.c_str(),
                                             map_ConfigToEffAndBkgd[{5000,*it_Config2}].first,
                                             map_ConfigToEffAndBkgd[{5000,*it_Config2}].second));
  }
  leg_EfficiencyVTimeWindow->Draw();
  c_Global->Print("Results.pdf");
  c_Global->Print("Results.pdf]");
  exit(0);
  int color = vec_Colors.at(globalIt % vec_Colors.size());
  h_Alex_FakeRateVNClusters   ->SetLineColor(color);
  h_Alex_FakeRateVNClustersLow->SetLineColor(color);
  h_Alex_EfficiencyVEvents    ->SetLineColor(color);
  h_Alex_EfficiencyVDistance  ->SetLineColor(color);
  h_Alex_EffGalaxy            ->SetLineColor(color);

  h_Alex_FakeRateVNClusters   ->SetLineStyle(2);
  h_Alex_FakeRateVNClustersLow->SetLineStyle(2);
  h_Alex_EfficiencyVEvents    ->SetLineStyle(2);
  h_Alex_EfficiencyVDistance  ->SetLineStyle(2);
  h_Alex_EffGalaxy            ->SetLineStyle(2);

  h_Alex_FakeRateVNClusters   ->SetLineWidth(2);
  h_Alex_FakeRateVNClustersLow->SetLineWidth(2);
  h_Alex_EfficiencyVEvents    ->SetLineWidth(2);
  h_Alex_EfficiencyVDistance  ->SetLineWidth(2);
  h_Alex_EffGalaxy            ->SetLineWidth(2);

  g_Alex_ROC                  ->SetMarkerColor(color);
  h_Alex_FakeRateVNClusters   ->SetMarkerColor(color);
  h_Alex_FakeRateVNClustersLow->SetMarkerColor(color);
  h_Alex_EfficiencyVEvents    ->SetMarkerColor(color);
  h_Alex_EfficiencyVDistance  ->SetMarkerColor(color);
  h_Alex_EffGalaxy            ->SetMarkerColor(color);
  g_Alex_ROC                  ->SetMarkerColor(color);

  
  


  THStack *stk_FakeRateVNClusters = new THStack("stk_FakeRateVNClusters", "Number of Clusters in Time Window Required to Trigger vs. Trigger Rate");
  TLegend *leg_FakeRateVNClusters = new TLegend(0.05, 0.05, 0.95, 0.95);
  leg_FakeRateVNClusters->SetTextSize(0.023);
  leg_FakeRateVNClusters->SetHeader(legHeader.c_str());

  for(it_Config  = map_ConfigToEffAndBkgd.begin();
      it_Config != map_ConfigToEffAndBkgd.end(); ++it_Config){
    stk_FakeRateVNClusters->Add(map_h_FakeRateVNClusters     [it_Config->first]);
    stk_FakeRateVNClusters->Add(map_h_FakeRateVNClustersLow  [it_Config->first]);
    leg_FakeRateVNClusters->AddEntry(map_h_FakeRateVNClusters[it_Config->first],
                                     Form(legEntryFormat.c_str(), it_Config->first.first, it_Config->second.first, it_Config->second.second), "L");
  }
  stk_FakeRateVNClusters->Add(h_Alex_FakeRateVNClusters);
  stk_FakeRateVNClusters->Add(h_Alex_FakeRateVNClustersLow);
  // leg_FakeRateVNClusters->AddEntry(h_Alex_FakeRateVNClusters, Form("Alex TW: 10000 ms, Eff: %.2f & Bkgd rate: %.2f Hz",
  //                                                                  map_ConfigToEffAndBkgd2.begin()->second.first,
  //                                                                  map_ConfigToEffAndBkgd2.begin()->second.second));

  c_Global->Draw();

  leg_FakeRateVNClusters->Draw();
  c_Global->Print("Results.pdf");

  c_Global->SetLogy();
  stk_FakeRateVNClusters->SetMinimum(1e-9);
  stk_FakeRateVNClusters->Draw("NOSTACK HIST");
  
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
  //leg_FakeRateVNClusters->Draw();
  c_Global->Print("Results.pdf");

  THStack *stk_EfficiencyVEvents = new THStack("stk_EfficiencyVEvents", "Efficiency vs. Number of Events in SN Burst, Fake Trigger Rate: 1/Month");
  TLegend *leg_EfficiencyVEvents = new TLegend(0.53, 0.65, 0.85, 0.85);
  leg_EfficiencyVEvents->SetTextSize(0.023);
  
  leg_EfficiencyVEvents->SetHeader(legHeader.c_str());
  for(it_Config  = map_ConfigToEffAndBkgd.begin();
      it_Config != map_ConfigToEffAndBkgd.end(); ++it_Config){
    stk_EfficiencyVEvents->Add(map_h_EfficiencyVEvents     [it_Config->first]);
    leg_EfficiencyVEvents->AddEntry(map_h_EfficiencyVEvents[it_Config->first], Form(legEntryFormat.c_str(), it_Config->first.first, it_Config->second.first, it_Config->second.second), "L");
  }
  stk_EfficiencyVEvents->Add(h_Alex_EfficiencyVEvents);
  leg_EfficiencyVEvents->AddEntry(h_Alex_EfficiencyVEvents, Form("Alex TW: 10000 ms, Eff: %.2f & Bkgd rate: %.2f Hz",
                                                                  map_ConfigToEffAndBkgd2.begin()->second.first,
                                                                  map_ConfigToEffAndBkgd2.begin()->second.second));
  c_Global->Clear();
  c_Global->SetLogy(false);
  c_Global->Draw();
  c_Global->SetLogx();
  stk_EfficiencyVEvents->Draw("NOSTACK");
  stk_EfficiencyVEvents->GetXaxis()->SetTitle("Number of Events in SN Burst");
  stk_EfficiencyVEvents->GetYaxis()->SetTitle("Efficiency");
  gPad->RedrawAxis();
  //leg_EfficiencyVEvents->Draw();
  c_Global->Print("Results.pdf");

  THStack *stk_EfficiencyVDistance = new THStack("stk_EfficiencyVDistance", "Efficiency vs. Distance to SN, Fake Trigger Rate: 1/Month");
  TLegend *leg_EfficiencyVDistance = new TLegend(0.15, 0.15, 0.48, 0.35);
  leg_EfficiencyVDistance->SetTextSize(0.023);

  leg_EfficiencyVDistance->SetHeader(legHeader.c_str());
  for(it_Config  = map_ConfigToEffAndBkgd.begin();
      it_Config != map_ConfigToEffAndBkgd.end(); ++it_Config){
    stk_EfficiencyVDistance->Add(map_h_EfficiencyVDistance[it_Config->first]);
    leg_EfficiencyVDistance->AddEntry(map_h_EfficiencyVDistance[it_Config->first], 
                                      Form(legEntryFormat.c_str(), it_Config->first.first, it_Config->second.first, it_Config->second.second), "L");
  }
  stk_EfficiencyVDistance->Add(h_Alex_EfficiencyVDistance);
  leg_EfficiencyVDistance->AddEntry(h_Alex_EfficiencyVDistance, Form("Alex TW: 10000 ms, Eff: %.2f & Bkgd rate: %.2f Hz",
                                                                     map_ConfigToEffAndBkgd2.begin()->second.first,
                                                                     map_ConfigToEffAndBkgd2.begin()->second.second));
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
  TLegend *leg_EffGalaxy = new TLegend(0.15, 0.15, 0.48, 0.35);
  leg_EffGalaxy->SetTextSize(0.023);

  leg_EffGalaxy->SetHeader(legHeader.c_str());
  stk_EffGalaxy->Add(h_SNProbabilityVDistance);

  for(it_Config =  map_ConfigToEffAndBkgd.begin();
      it_Config != map_ConfigToEffAndBkgd.end(); ++it_Config){
    stk_EffGalaxy->Add(map_h_EffGalaxy[it_Config->first]);
    leg_EffGalaxy->AddEntry(map_h_EffGalaxy[it_Config->first],
                            Form(legEntryFormat.c_str(), it_Config->first.first, it_Config->second.first, it_Config->second.second), "L");
  }
  stk_EffGalaxy->Add(h_Alex_EffGalaxy);
  leg_EffGalaxy->AddEntry(h_Alex_EffGalaxy, Form("Alex TW: 10000 ms, Eff: %.2f & Bkgd rate: %.2f Hz",
                                                 map_ConfigToEffAndBkgd2.begin()->second.first,
                                                 map_ConfigToEffAndBkgd2.begin()->second.second));
  leg_EffGalaxy->AddEntry(h_SNProbabilityVDistance, "SN Probability", "L");
  c_Global->Clear();
  c_Global->Draw();
  c_Global->SetLogx(false);
  c_Global->SetLogy();
  stk_EffGalaxy->Draw("NOSTACK");
  stk_EffGalaxy->GetXaxis()->SetTitle("SN Distance, (kpc)");
  stk_EffGalaxy->GetYaxis()->SetTitle("Efficiency x SN Probability");
  stk_EffGalaxy->Draw("NOSTACK");
  gPad->RedrawAxis();
  //leg_EffGalaxy->Draw();
  c_Global->Print("Results.pdf");


  TLegend *leg_ROC    = new TLegend(0.15, 0.68, 0.48, 0.88);
  leg_ROC->SetTextSize(0.023);
  leg_ROC->SetHeader(legHeader.c_str());
  c_Global->Clear();
  c_Global->Draw();
  double minX=0.9, maxX=1;
  double minY=10e-15, maxY=20e2;

  map_g_ROC.begin()->second->GetXaxis()->SetLimits(minX, maxX);
  map_g_ROC.begin()->second->GetYaxis()->SetLimits(minY, maxY);
  map_g_ROC.begin()->second->SetTitle("Fake Trigger Rate vs. Galactic Neighbourhood Coverage");
  map_g_ROC.begin()->second->GetXaxis()->SetTitle("Galactic Neighbourhood Coverage");
  map_g_ROC.begin()->second->GetYaxis()->SetTitle("Fake Trigger Rate, (Hz)");
  map_g_ROC.begin()->second->Draw("AP");

  for(it_Config =  map_ConfigToEffAndBkgd.begin();
      it_Config != map_ConfigToEffAndBkgd.end(); ++it_Config){
    leg_ROC->AddEntry(map_g_ROC[it_Config->first], Form(legEntryFormat.c_str(), it_Config->first.first, it_Config->second.first, it_Config->second.second), "P");
    map_g_ROC[it_Config->first]->Draw("P");
  }
  g_Alex_ROC->Draw("P");
  leg_ROC->AddEntry(g_Alex_ROC, Form("Alex TW: 10000 ms, Eff: %.2f & Bkgd rate: %.2f Hz",
                                     map_ConfigToEffAndBkgd2.begin()->second.first,
                                     map_ConfigToEffAndBkgd2.begin()->second.second));
  TLine *l_perMonth_2 = new TLine(minX, 4.13e-7, maxX, 4.13e-7);
  l_perMonth_2->SetLineColor(1);
  l_perMonth_2->SetLineWidth(3);
  TText *t_perMonth_2 = new TText(minX+0.015, 8e-7, "1/Month");
  gPad->RedrawAxis();
  //leg_ROC->Draw();
  l_perMonth_2->Draw();
  t_perMonth_2->Draw();
  c_Global->Print("Results.pdf");
  c_Global->Print("Results.pdf]");

  return 0;
}
