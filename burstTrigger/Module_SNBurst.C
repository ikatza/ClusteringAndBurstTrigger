#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <TString.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TMath.h>
#include <TGraph.h>

using namespace std;

TFile *f_Output;
TFile *f_Theory;
int threshold = 10;
const bool reproduceAlexResult = false;
const bool faster = true;

std::map<double,int> makeDistanceToEventsNeighbourhoodMap(TH1D* const &h_SNProbabilityVDistance,
                                                          TF1* const &f_EventsVSNDistance_10kt){
  std::map<double,int> map_DistanceToEvents_Neighbourhood;
  for(int i = 1; i < h_SNProbabilityVDistance->GetSize()-1; i++){
    double distance = h_SNProbabilityVDistance->GetBinCenter(i);
    int    nEvents  = std::floor(f_EventsVSNDistance_10kt->Eval(distance));
    map_DistanceToEvents_Neighbourhood[distance] = nEvents;
  }

  return map_DistanceToEvents_Neighbourhood;
}


Double_t PoissonIntegral(Double_t x, Double_t mean){
  double value = 0;
  for (int i=0; i<x; ++i) {
    value += TMath::Power(mean,i) / TMath::Factorial(i);
  }
  value *=TMath::Exp(-mean);
  return value;
}


void makeFakeRateVNClusters(double const &mean, double const &rms, double const &bkgd,
                            int const &config,double const &timeW,
                            std::map<double,double> &map_FakeRateToNClusters){

  int    nBins  = 2000;
  double rmsMax = 10*rms;
  TH1D  *h_FakeRateVNClusters    = new TH1D(Form("h_FakeRateVNClusters_Config%i_TW%.3f",config, timeW),
                                            Form("h_FakeRateVNClusters_Config%i_TW%.3f",config, timeW), nBins,  mean, mean+rmsMax);
  TH1D  *h_FakeRateVNClustersLow = new TH1D(Form("h_FakeRateVNClustersLow_Config%i_TW%.3f",config, timeW),
                                            Form("h_FakeRateVNClustersLow_Config%i_TW%.3f",config, timeW), nBins, 0, mean);

  for(int i = 1; i < h_FakeRateVNClusters->GetSize()-1; i++){
    double nClusters = h_FakeRateVNClusters->GetBinCenter(i);
    if(nClusters<0){
      std::cout << "Clusters < 0 " << std::endl;
      continue;
    }else if(mean > threshold){
      double nRMS      = (nClusters-mean)/rms;
      double fraction  = (1-TMath::Erf(nRMS/std::sqrt(2.)))/2.;
      double fakeRate  = fraction*bkgd;

      h_FakeRateVNClusters->SetBinContent(i, fakeRate);
      map_FakeRateToNClusters[fakeRate] = nClusters;
    }else{
      double integral = 1-PoissonIntegral(nClusters,mean);
      double fakeRate  = integral*bkgd;

      h_FakeRateVNClusters->SetBinContent(i, fakeRate);
      map_FakeRateToNClusters[fakeRate] = nClusters;
    }
  }

  //GO BELOW THE MEAN.
  for(int i = 1; i < h_FakeRateVNClustersLow->GetSize()-1; i++){
    double nClusters = h_FakeRateVNClustersLow->GetBinCenter(i);
    if(nClusters<0){
      continue;
    }else if(mean > threshold){
      double nRMS      = (nClusters-mean)/rms;
      double fraction  = (1-TMath::Erf(nRMS/std::sqrt(2.)))/2.;
      double fakeRate  = (fraction)*bkgd;

      h_FakeRateVNClustersLow->SetBinContent(i, fakeRate);
    }else{
      double integral = 1-PoissonIntegral(nClusters,mean);
      double fakeRate  = integral*bkgd;

      h_FakeRateVNClustersLow->SetBinContent(i, fakeRate);
    }
  }
  f_Output->cd();
  h_FakeRateVNClusters->Write();
  h_FakeRateVNClustersLow->Write();

  delete h_FakeRateVNClusters;
  delete h_FakeRateVNClustersLow;
  h_FakeRateVNClusters = NULL;
  h_FakeRateVNClustersLow = NULL;

  return;
}


std::map<std::pair<double,int>,double>  makeEfficiencyVEvents(double const &eff,       double const &mean,     double const &rms, double const &fracInTW,
                                                              double const &nClusters, int    const &burstMin, int    const &burstMax,
                                                              std::map<double,int> &map_ClustersToMaxEffEvent){
  std::map<std::pair<double,int>,double> map_NClustersAndEventsToBurstEff;
  bool isOne = false;
  double f = fracInTW;

  for(int i = burstMin; i <= burstMax; i++){
    double integral      = 0;
    int    totalClusters = (int)std::ceil(mean + double(i)*f);
    TF1 *f_Poisson = new TF1("f_Poisson", "TMath::PoissonI(x,[0])", 0, mean + 20*rms);
    int  nptx = std::ceil(TMath::Abs(10*(mean + 20*rms)));
    f_Poisson->SetNpx(nptx>100000?100000:nptx);
    f_Poisson->SetParameter(0,totalClusters);

    if(isOne == false)
    {
      integral = 1. - f_Poisson->Integral(0, nClusters/eff, 1e-4);

      if(integral>=0 && integral<=1)
      {
        map_NClustersAndEventsToBurstEff[{nClusters,i}] = integral;
        if(integral==1)
        {
          map_ClustersToMaxEffEvent[nClusters] = i;
          break;
        }
      }
      else if(integral<0)
      {
        map_NClustersAndEventsToBurstEff[{nClusters,i}] = 0;
      }
    }

    if(f_Poisson){
      delete f_Poisson;
      f_Poisson = NULL;
    }
  }

  return map_NClustersAndEventsToBurstEff;
}


void makeNeighbourhoodEfficiency(TH1D* const &h_SNProbabilityVDistance, std::map<double,int> &map_DistanceToEvents_Neighbourhood,
                                 std::map<std::pair<double,int>,double> &map_NClustersAndEventsToBurstEff, double const &nClusters,
                                 TH1D* &h_NeighbourhoodEffiency, std::map<double,int> &map_ClustersToMaxEffEvent,
                                 std::map<double,double> &map_ClustersToCoverage)
{
  for(int i = 1; i < h_SNProbabilityVDistance->GetSize()-1; i++)
  {
    double prob     = h_SNProbabilityVDistance->GetBinContent(i);
    double distance = h_SNProbabilityVDistance->GetBinCenter(i);
    int    nEvents  = map_DistanceToEvents_Neighbourhood[distance];
    double probXEff, burstEff;
    if(nEvents>=map_ClustersToMaxEffEvent[nClusters])
    {
      probXEff = prob*1;
    }
    else
    {
      burstEff = map_NClustersAndEventsToBurstEff[{nClusters,nEvents}];
      probXEff = prob*burstEff;
    }
    h_NeighbourhoodEffiency->SetBinContent(i,probXEff);
  }

  double integral = h_NeighbourhoodEffiency->Integral();
  map_ClustersToCoverage[nClusters] = integral;

  return;
}




int main()
{
  //GRAB INFORMATION FROM THE CLUSTERING ALGORITHM AND DEFINE OUTPUT FILES.
  TString s_Filename = "GH_SNMC";
  ifstream inFile;
  inFile.open("Analyse_"+s_Filename+".txt");
  ofstream outFile;
  outFile.open("Analyse_SNBurst_"+s_Filename+".txt");

  //DEFINE PARAMETERS.
  int burstMin(1), burstMax(30e4);

  double cut_PerYear = 3.44e-8;

  //GRAB THEORY PLOTS.
  f_Theory = new TFile("SNTheoryDistributions.root","READ");
  f_Output = new TFile("Analyse_SNBurst_"+s_Filename+".root", "RECREATE");
  TFile *f_TimeProfile = new TFile("TimeProfile.root", "READ");
  TH1D* h_TimeProfile = (TH1D*)f_TimeProfile->Get("h_MarlTime");

  TH1D *h_SNProbabilityVDistance = (TH1D*)f_Theory->Get("h_SNProbabilityVDistance_LMC");
  h_SNProbabilityVDistance->SetDirectory(0);

  TF1  *f_EventsVSNDistance_10kt = (TF1*)f_Theory->Get("f_EventsVSNDistance_10kt_100kpc");
  double gradient  = f_EventsVSNDistance_10kt->GetParameter(0);
  double intercept = f_EventsVSNDistance_10kt->GetParameter(1);
  TF1 *f_Inverse   = new TF1("f_Inverse", "TMath::Power(x/(TMath::Power(10,[0])),1/[1])", 1,40e4);
  f_Inverse->SetParameter(0, intercept);
  f_Inverse->SetParameter(1, gradient);
  double min_Distance = f_Inverse->Eval(burstMax);
  double max_Distance = f_Inverse->Eval(burstMin);
  std::map<double,int> map_DistanceToEvents_Neighbourhood = makeDistanceToEventsNeighbourhoodMap(h_SNProbabilityVDistance,
                                                                                                 f_EventsVSNDistance_10kt);

  std::map<std::pair<int,double>,std::pair<double,double>> map_ConfigToEffAndBkgd;
  int Config;
  double TimeWindow, Eff, Bkgd;

  while(inFile >> Config >> TimeWindow >> Eff >> Bkgd){
    map_ConfigToEffAndBkgd[std::make_pair(Config,TimeWindow)] = {Eff,Bkgd};
  }

  //LOOP AROUND THE CLUSTERING CONFIGURATIONS.
  std::map<int,std::pair<double,double>>::iterator it_ConfigToEffAndBkgd;
  for(auto const& it_ConfigToEffAndBkgd : map_ConfigToEffAndBkgd){

    int    Config = it_ConfigToEffAndBkgd.first.first;
    double TimeW  = it_ConfigToEffAndBkgd.first.second;
    double eff    = it_ConfigToEffAndBkgd.second.first;
    double bkgd   = it_ConfigToEffAndBkgd.second.second;
    //ASSUME THE NUMBER OF CLUSTERS IN A GIVEN TIME WINDOW IS GAUSSIAN ABOUT THE GIVEN BACKGROUND RATE WITH AN RMS OF SQRT(MEAN).
    double mean = bkgd*TimeW;
    double rms  = std::sqrt(mean);
    double fracInTW = h_TimeProfile->Integral(0,h_TimeProfile->FindBin(TimeW*1000));
    if (reproduceAlexResult) {
      threshold = 0;
      fracInTW = 1;
    }
    std::cout << "Frac " << fracInTW << std::endl;
    if(fracInTW<0 || fracInTW>1){
      std::cout << "fracInTWh is "  << fracInTW << " !!!" << std::endl;
      exit(1);
    }
    std::cout << " sig/mean " << eff*fracInTW/mean << std::endl;
    std::map<double,double> map_FakeRateToNClusters;
    makeFakeRateVNClusters(mean, rms, bkgd, Config,TimeW, map_FakeRateToNClusters);

    //PULL OUT THE PER YEAR CLUSTER CUT FOR THIS CONFIGURATION.
    double cut_KeyNClustersForOnePerYear = 0.;
    double cut_NClustersForOnePerYear    = 0.;
    double lastKey=-1, lastN=-1;
    for(auto const& it_FakeRateToNClusters: map_FakeRateToNClusters){
      if(it_FakeRateToNClusters.first >= cut_PerYear){
        cut_KeyNClustersForOnePerYear = lastKey;
        cut_NClustersForOnePerYear    = lastN;
        break;
      }
      lastKey = it_FakeRateToNClusters.first;
      lastN   = it_FakeRateToNClusters.second;
    }
    if(cut_KeyNClustersForOnePerYear == -1){
      std::cout << "KeyNClustersForOnePerYear is " << cut_KeyNClustersForOnePerYear << std::endl;
    } else {
      outFile   << Config << " "
                << TimeW << " "
                << eff << " "
                << bkgd << " "
                << cut_NClustersForOnePerYear << std::endl;
    }
    std::cout << "CONFIG: " << Config << " TIMEWINDOW: " << TimeW
              << ", EFF: " << eff
              << ", BKGD: " << Bkgd
              << ", PERYEAR CUT:  " << cut_NClustersForOnePerYear << std::endl;

    //LOOP OVER THE DIFFERENT FAKE RATES AND CUTS TO GET FAKE RATE
    //AGAINST GALACTIC COVERAGE.
    std::map<std::pair<double,int>,double> map_NClustersAndEventsToBurstEff;
    std::map<double,double> map_ClustersToCoverage;
    std::map<double,int>    map_ClustersToMaxEffEvent;
    int count_Loop = 0;
    // THE fraction of events that are expected to fall in the time window
    std::cout << "BIN: " << h_TimeProfile->FindBin((double)TimeW) << std::endl;
    std::cout << "map_FakeRateToNClusters.size() " << map_FakeRateToNClusters.size() << std::endl;
    for(auto const& it_FakeRateToNClusters : map_FakeRateToNClusters){
      if(count_Loop % 500 == 0)
        std::cout << "WORKING ON CONFIG: " << Config << " TIMEWINDOW: " << TimeW
                  << ", ITERATION: " << count_Loop << std::endl;

      TH1D *h_NeighbourhoodEffiency = (TH1D*)h_SNProbabilityVDistance->Clone();
      //h_NeighbourhoodEffiency->Reset();
      h_NeighbourhoodEffiency->SetName(Form("h_NeighbourhoodEffiency_Config%i_TW%0.3f",Config,TimeW));
      h_NeighbourhoodEffiency->SetDirectory(0);
      // FOR EACH BURST SIZE
      if(!faster ||
         (faster && it_FakeRateToNClusters.first  == cut_KeyNClustersForOnePerYear && it_FakeRateToNClusters.second == cut_NClustersForOnePerYear)) {
        map_NClustersAndEventsToBurstEff = makeEfficiencyVEvents(eff, mean, rms, fracInTW, it_FakeRateToNClusters.second,
                                                                 burstMin, burstMax,
                                                                 map_ClustersToMaxEffEvent);
        makeNeighbourhoodEfficiency(h_SNProbabilityVDistance, map_DistanceToEvents_Neighbourhood,
                                    map_NClustersAndEventsToBurstEff, it_FakeRateToNClusters.second,
                                    h_NeighbourhoodEffiency, map_ClustersToMaxEffEvent, map_ClustersToCoverage);
      }

      // MAKE THE 1 PER YEAR HISTOGRAMS.
      if(it_FakeRateToNClusters.first  == cut_KeyNClustersForOnePerYear &&
         it_FakeRateToNClusters.second == cut_NClustersForOnePerYear){
        std::cout << "Iteration " << count_Loop << std::endl;
        std::cout << "it_FakeRateToNClusters->first  " << it_FakeRateToNClusters.first  << " cut_KeyNClustersForOnePerYear " << cut_KeyNClustersForOnePerYear << std::endl;
        std::cout << "it_FakeRateToNClusters->second " << it_FakeRateToNClusters.second << " cut_NClustersForOnePerYear    " << cut_NClustersForOnePerYear    << std::endl;
        TH1D *h_EfficiencyVEvents   = new TH1D(Form("h_EfficiencyVEvents_Config%i_TW%.3f",Config,TimeW),
                                               Form("h_EfficiencyVEvents_Config%i_TW%.3f",Config,TimeW),
                                               burstMax-burstMin+1, burstMin, burstMax);
        for(int i = burstMin; i <= burstMax; i++){
          int bin = h_EfficiencyVEvents->FindBin(i);
          h_EfficiencyVEvents->SetBinContent(bin,map_NClustersAndEventsToBurstEff[{cut_NClustersForOnePerYear,i}]);
        }
        f_Output->cd();
        h_EfficiencyVEvents->Write();

        TH1D *h_EfficiencyVDistance = new TH1D(Form("h_EfficiencyVDistance_Config%i_TW%.3f",Config,TimeW),
                                               Form("h_EfficiencyVDistance_Config%i_TW%.3f",Config,TimeW),
                                               burstMax-burstMin+1, min_Distance, max_Distance);
        for(int i = 1; i < h_EfficiencyVDistance->GetSize()-1; i++){
          double distance   = h_EfficiencyVDistance->GetBinCenter(i);
          double nEvents    = f_EventsVSNDistance_10kt->Eval(distance);
          int    nEventsBin = h_EfficiencyVEvents->FindBin(nEvents);
          double burstEff   = h_EfficiencyVEvents->GetBinContent(nEventsBin);

          h_EfficiencyVDistance->SetBinContent(i,burstEff);
        }
        f_Output->cd();
        h_EfficiencyVDistance->Write();
        if(h_EfficiencyVEvents){
          delete h_EfficiencyVEvents;
          h_EfficiencyVEvents = NULL;
        }
        if(h_EfficiencyVDistance){
          delete h_EfficiencyVDistance;
          h_EfficiencyVDistance = NULL;
        }
        h_NeighbourhoodEffiency->Write();
      }

      if(h_NeighbourhoodEffiency){
        delete h_NeighbourhoodEffiency;
        h_NeighbourhoodEffiency = NULL;
      }

      count_Loop++;
    }

    //EXTRACT THE INFORMATION WE HAVE PICKED UP AND MAKE THE ROC PLOT.
    TGraph *g_ROC = new TGraph(map_FakeRateToNClusters.size());
    g_ROC->SetName(Form("g_ROC_Config%i_TW%0.3f",Config,TimeW));
    count_Loop = 0;

    for(auto const& it_FakeRateToNClusters : map_FakeRateToNClusters){
      //if(map_ClustersToCoverage.count(it_FakeRateToNClusters.second) == 1){
        // std::cout << "it_FakeRateToNClusters.first " << it_FakeRateToNClusters.first << std::endl;
        // std::cout << "it_FakeRateToNClusters.second " << it_FakeRateToNClusters.second << std::endl;
        // std::cout << "map_ClustersToCoverage[it_FakeRateToNClusters.second] " << map_ClustersToCoverage[it_FakeRateToNClusters.second] << std::endl;
        g_ROC->SetPoint(count_Loop,map_ClustersToCoverage[it_FakeRateToNClusters.second], it_FakeRateToNClusters.first);
        count_Loop++;
        //}
    }
    f_Output->cd();
    g_ROC->Write();
  }
  f_Output->Close();
  f_Theory->Close();
  return 0;
}
