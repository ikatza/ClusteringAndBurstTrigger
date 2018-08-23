#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <iomanip>

#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TString.h"
#include "TROOT.h"

#include "Clustering.hh"
#include <unistd.h>

std::map<int,double> map_TypeToWeight;

double getWeight(const int ClusterType) {
  return map_TypeToWeight[ClusterType];
};


int getClusterType(const std::vector<int>* vec_GenType) {
  int type = -1;

  std::map<int,int> map_GenTypeToCount;
  for(auto const & it : (*vec_GenType))
    map_GenTypeToCount[it]++;  

  int largestCount = 0;
  for(auto const& it : map_GenTypeToCount) {
    if(it.second>largestCount && it.first!=1) {
      largestCount = it.second;
      type         = it.first; 
    }
  }

  return type;
};


void FillEventCountMap(TTree* ClusteredHit,
                       size_t nHitCut,
                       std::map<int,double>& map_Config_nBackCluster,
                       std::map<int,std::map<int,std::map<int,int>>>& map_ConfigEventIndex_nSignCluster){
  
  int               in_Event;
  int               in_Config;
  double            in_Type;
  int               in_MarleyIndex;
  std::vector<int>* in_GenType = NULL;
  
  ClusteredHit->SetBranchAddress("Event",        &in_Event      );
  ClusteredHit->SetBranchAddress("Config",       &in_Config     );
  ClusteredHit->SetBranchAddress("Type",         &in_Type       );
  ClusteredHit->SetBranchAddress("MarleyIndex",  &in_MarleyIndex);
  ClusteredHit->SetBranchAddress("MarleyIndex",  &in_MarleyIndex);
  if (ClusteredHit->GetListOfBranches()->FindObject("GenType")) {
    ClusteredHit->SetBranchAddress("GenType", &in_GenType);
  } else {
    ClusteredHit->SetBranchAddress("Hit_GenType", &in_GenType);
  }
  for(int i = 0; i < ClusteredHit->GetEntries(); i++) {
    ClusteredHit->GetEntry(i);
    if (in_GenType->size() <= nHitCut) continue;
    if (in_Type == 1) {
      map_ConfigEventIndex_nSignCluster[in_Config][in_Event][in_MarleyIndex]++;
    } else {
      int    ClusterType = getClusterType(in_GenType);
      double weight = getWeight(ClusterType);
      map_Config_nBackCluster[in_Config] += weight;
    }
  }
  return;
  
};


void GetEfficiency(const int Config,
                   std::map<int,std::map<int,std::map<int,int>>> const& map_ConfigEventIndex_nSignCluster,
                   std::map<int,int> const& nMarleyEvent,
                   std::pair<double,double>& efficiency) {

  for(auto const& it0 : map_ConfigEventIndex_nSignCluster) {
    int ThisConfig = it0.first;
    if (ThisConfig != Config) continue;
    std::map<int, std::map<int,int>> const& map_EventIndex_nSignCluster = it0.second;
    efficiency.first = 0;
    for(auto const& it1 : map_EventIndex_nSignCluster) {
      int ThisEvent = it1.first;
      std::map<int,int> const& map_Index_nSignCluster = it1.second;
      for (int iMarley=0; iMarley<nMarleyEvent.at(ThisEvent); ++iMarley) {
        if (map_Index_nSignCluster.find(iMarley) != map_Index_nSignCluster.end())
          if (map_Index_nSignCluster.at(iMarley) > 0)
            efficiency.first++;
      }
    }
    double nMarley=0;
    for(auto const& itEvent : nMarleyEvent)
      nMarley += itEvent.second;
    std::cout << efficiency.first << " signal events in config " << Config
              << " (" << nMarley << " were generated)." << std::endl;
    efficiency.second = sqrt(1./efficiency.first+1./((double)nMarley));
    efficiency.first /=(double)nMarley;
    efficiency.second *= efficiency.first; 
    return;
  }
  
};

void GetBackgroundRate(const int Config,
                       const std::map<int,double> map_Config_nBackCluster,
                       const int nEvent,
                       const double EventTime,
                       const double DetectorScaling,
                       std::pair<double,double>& BackgroundRate) {

  BackgroundRate.first = map_Config_nBackCluster.at(Config) / ((double) nEvent * EventTime) / (DetectorScaling); 
  BackgroundRate.second = sqrt(1./map_Config_nBackCluster.at(Config)+1./((double)nEvent)) * BackgroundRate.first;
  return;
};

int main(int argc, char** argv) {
  int opt;
  extern char *optarg;
  extern int   optopt;

  std::string InputFile  = "";
  std::string OutputFile = "";
  int Config = -1;
  int nHitCut = 0;
  double DetectorScaling = 0.12;
  while ( (opt = getopt(argc, argv, "i:o:c:s:h:")) != -1 ) {  // for each option...
    switch ( opt ) {
    case 'i':
      InputFile  = std::string(optarg);
      break;
    case 'o':
      OutputFile = std::string(optarg);
      break;
    case 'c':
      Config = atoi(optarg);
      break;
    case 's':
      DetectorScaling = atoi(optarg);
      break;
    case 'h':
      std::cout << optarg<< std::endl;
      nHitCut = atoi(optarg);
      break;
    case '?':  // unknown option...
      std::cerr << "Unknown option: '" << char(optopt) << "'!" << std::endl;
      break;
    }
  }

  Clustering c;
  if (Config >= c.GetNConfig()) {
    std::cerr << "Requested config (" << Config << ") doesn't exit (Config can be from 0 to " << c.GetNConfig()-1 << std::endl;
    exit(1);
  }

  int nConfig = (Config==-1? c.GetNConfig(): Config+1);
  int iterConfig = (Config==-1? 0:Config);

  TFile *f_Input = new TFile(InputFile.c_str(), "READ");
  TTree *ClusteredWireHit    = (TTree*)f_Input->Get("ClusteredWireHit");
  TTree *ClusteredOpticalHit = (TTree*)f_Input->Get("ClusteredOpticalHit");
  TTree *TrueInfo            = (TTree*)f_Input->Get("TrueInfo");
  
  int nEventOriginally = (int)TrueInfo->GetEntries();
  std::cout << "THERE WERE " << nEventOriginally << " EVENTS IN THE ORIGINAL SAMPLE AND " << nConfig << " CONFIGURATIONS" << std::endl;
 
  map_TypeToWeight[0] = 1.; //NOISE 
  map_TypeToWeight[1] = 1.; //MARLEY
  map_TypeToWeight[2] = 1.; //CO60
  map_TypeToWeight[3] = 1.; //K40
  map_TypeToWeight[4] = 1.; //AR39
  map_TypeToWeight[5] = 1.; //n
  map_TypeToWeight[6] = 1.; //KR
  map_TypeToWeight[7] = 1.; //PO
  map_TypeToWeight[8] = 1.; //RN
  map_TypeToWeight[9] = 1.; //AR42

  std::vector<double>* dummy=NULL;

  int nEvent = TrueInfo->GetEntries();
  int in_Event;
  TrueInfo->SetBranchAddress("Event", &in_Event);
  TrueInfo->SetBranchAddress("MarlTime", &dummy);

  std::map<int,int> map_Event_nMarley;
  int nMarleyEvent=0;
  for(int i = 0; i < nEvent; i++) {
    TrueInfo->GetEntry(i);
    map_Event_nMarley[in_Event] = (int)dummy->size();
    nMarleyEvent += (int)dummy->size();
  }
  std::cout << "There were " << nMarleyEvent << " Marley events in " << map_Event_nMarley.size() << " Larsoft events." << std::endl;
  //OVERALL EFFICIENCIES AND BACKGROUND RATES.
  std::map<int,double> map_Config_nWireBackCluster;
  std::map<int,std::map<int,std::map<int,int>>> map_ConfigEventIndex_nWireSignCluster;
  FillEventCountMap(ClusteredWireHit,
                    nHitCut,
                    map_Config_nWireBackCluster,
                    map_ConfigEventIndex_nWireSignCluster);

  std::map<int,double> map_Config_nOpticalBackCluster;
  std::map<int,std::map<int,std::map<int,int>>> map_ConfigEventIndex_nOpticalSignCluster;
  FillEventCountMap(ClusteredOpticalHit,
                    nHitCut,
                    map_Config_nOpticalBackCluster,
                    map_ConfigEventIndex_nOpticalSignCluster);

  std::map<int,std::pair<double,double>> map_Config_WireEff;
  std::map<int,std::pair<double,double>> map_Config_WireBackRate;
  std::map<int,std::pair<double,double>> map_Config_OpticalEff;
  std::map<int,std::pair<double,double>> map_Config_OpticalBackRate;

  for (int c=iterConfig; c<nConfig; ++c) {
    std::cout << "N background wire clusters " << map_Config_nWireBackCluster[c] << " in config " << c << std::endl;
    std::cout << "N background opti clusters " << map_Config_nOpticalBackCluster[c] << " in config " << c << std::endl;
    GetEfficiency(c, map_ConfigEventIndex_nWireSignCluster, map_Event_nMarley,
                  map_Config_WireEff[c]);

    GetEfficiency(c, map_ConfigEventIndex_nOpticalSignCluster, map_Event_nMarley,
                  map_Config_OpticalEff[c]);

    GetBackgroundRate(c, map_Config_nWireBackCluster,nEventOriginally, 2.246e-3, DetectorScaling,
                      map_Config_WireBackRate[c]);

    GetBackgroundRate(c, map_Config_nOpticalBackCluster,nEventOriginally, 3.2e-3*3., DetectorScaling,
                      map_Config_OpticalBackRate[c]);
  }
  
  std::ofstream txt_Result;
  txt_Result.open(OutputFile.c_str());
  for(auto const& it : map_Config_WireEff) {
    int c = it.first;
    std::cout << std::setprecision(4);
    std::cout  << "Config (Wire) " << c
               << ", SN efficiency: " << it.second.first << " +/- " << it.second.second
               << " Background rate in 10kt (Hz): " << map_Config_WireBackRate[c].first << " +/- " << map_Config_WireBackRate[c].second << std::endl;
    txt_Result << c
               << ", " << it.second.first << ", " << it.second.second
               << ", " << map_Config_WireBackRate[c].first
               << ", " << map_Config_WireBackRate[c].second << std::endl; 
  }
  
  for(auto const& it : map_Config_WireEff) {
    int c = it.first;
    std::cout  << "Config (Opti) " << c
               << ", SN efficiency: " << map_Config_OpticalEff[c].first << " +/- " << map_Config_OpticalEff[c].second
               << " Background rate in 10kt (Hz): " << map_Config_OpticalBackRate[c].first << " +/- " << map_Config_OpticalBackRate[c].second << std::endl;
  }
  return 0;
}
