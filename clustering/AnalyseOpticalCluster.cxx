#include <algorithm>    // std::min_element, std::max_element
#include <iostream>
#include <fstream>
#include <unistd.h>


#include "Helper.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TPad.h"
#include "TTree.h"
#include "TLegend.h"
#include "TEfficiency.h"

TProfile* SetUpTProfileGenType(std::string name, std::string title) {
  Helper h;
  TProfile* p_ = new TProfile(name.c_str(), title.c_str(), h.GenName.size(), -0.5, (double)h.GenName.size()-0.5);
  for(auto const& it : h.ShortGenName)
    p_->GetXaxis()->SetBinLabel(it.first+1, it.second.c_str());
  return p_;

};

std::map<int, int> GetMapOfHit(std::vector<int> const* GenType){
  std::map<int,int> map_;
  for(std::vector<int>::const_iterator it=GenType->begin(); it!=GenType->end(); ++it)
    map_[*it]++;
  return map_;
};


int main(int argc, char** argv){

  int opt;
  // Shut GetOpt error messages down (return '?'):
  extern char *optarg;
  extern int  optopt;

  int nHitCut=-1;
  std::string InputFile="";
  std::string OutputFile="OpticalCluster.pdf";
  int InputConfig=5;
  while ( (opt = getopt(argc, argv, "n:o:i:c:")) != -1 ) {  // for each option...
    switch ( opt ) {
    case 'n':
      nHitCut = atoi(optarg);
      break;
    case 'i':
      InputFile = optarg;
      break;
    case 'o':
      OutputFile = optarg;
      break;
    case 'c':
        InputConfig = atoi(optarg);
        break;
    case '?':  // unknown option...
      std::cerr << "Unknown option: '" << char(optopt) << "'!" << std::endl;
      break;
    }
  }

  if (InputFile == "") {
    std::cout << "You need to provide an input file -i ClusteringOutput.root" << std::endl;
  }

  std::cout << "Saving output in " << OutputFile << std::endl;
  TFile *f_Input = new TFile(InputFile.c_str(), "READ");

  int    Cluster       ;
  int    Event         ;
  int    MarleyIndex   ;
  int    EventTrue     ;
  int    Config        ;
  int    NChan         ;
  int    Type          ;
  int    NHits         ;
  float  FirstTimeHit  ;
  float  LastTimeHit   ;
  float  TimeWidth     ;
  double RecClusterPosX;
  double RecClusterPosY;
  double RecClusterPosZ;
  double YWidth        ;
  double ZWidth        ;
  double SumPE         ;


  std::vector<int>    * PDSHit_GenType   = NULL;
  std::vector<double> * PDSHit_X         = NULL;
  std::vector<double> * PDSHit_Y         = NULL;
  std::vector<double> * PDSHit_Z         = NULL;
  std::vector<double> * PDSHit_PeakTime  = NULL;
  std::vector<double> * PDSHit_Width     = NULL;
  std::vector<double> * PDSHit_PE        = NULL;
  std::vector<double> * PDSHit_OpChannel = NULL;

  std::vector<double> * MarlTime = NULL;
  std::vector<double> * ENu      = NULL;
  std::vector<double> * ENu_Lep  = NULL;
  std::vector<double> * PosX     = NULL;
  std::vector<double> * PosY     = NULL;
  std::vector<double> * PosZ     = NULL;
  std::vector<double> * PosT     = NULL;
  std::vector<double> * DirX     = NULL;
  std::vector<double> * DirY     = NULL;
  std::vector<double> * DirZ     = NULL;

  TTree* t_ClusteredOptHit = (TTree*)f_Input->Get("ClusteredOpticalHit");
  t_ClusteredOptHit->SetBranchAddress("Cluster",        &Cluster       );
  t_ClusteredOptHit->SetBranchAddress("Config",         &Config        );
  t_ClusteredOptHit->SetBranchAddress("Event",          &Event         );
  t_ClusteredOptHit->SetBranchAddress("MarleyIndex",    &MarleyIndex   );
  t_ClusteredOptHit->SetBranchAddress("YWidth",         &YWidth        );
  t_ClusteredOptHit->SetBranchAddress("ZWidth",         &ZWidth        );
  t_ClusteredOptHit->SetBranchAddress("NChan",          &NChan         );
  t_ClusteredOptHit->SetBranchAddress("Type",           &Type          );
  t_ClusteredOptHit->SetBranchAddress("NHits",          &NHits         );
  t_ClusteredOptHit->SetBranchAddress("SumPE",          &SumPE         );
  t_ClusteredOptHit->SetBranchAddress("FirstTimeHit",   &FirstTimeHit  );
  t_ClusteredOptHit->SetBranchAddress("LastTimeHit",    &LastTimeHit   );
  t_ClusteredOptHit->SetBranchAddress("TimeWidth",      &TimeWidth     );
  t_ClusteredOptHit->SetBranchAddress("RecClusterPosX", &RecClusterPosX);
  t_ClusteredOptHit->SetBranchAddress("RecClusterPosY", &RecClusterPosY);
  t_ClusteredOptHit->SetBranchAddress("RecClusterPosZ", &RecClusterPosZ);
  t_ClusteredOptHit->SetBranchAddress("Hit_GenType"   , &PDSHit_GenType  );
  t_ClusteredOptHit->SetBranchAddress("Hit_X"         , &PDSHit_X        );
  t_ClusteredOptHit->SetBranchAddress("Hit_Y"         , &PDSHit_Y        );
  t_ClusteredOptHit->SetBranchAddress("Hit_Z"         , &PDSHit_Z        );
  t_ClusteredOptHit->SetBranchAddress("Hit_PeakTime"  , &PDSHit_PeakTime );
  t_ClusteredOptHit->SetBranchAddress("Hit_Width"     , &PDSHit_Width    );
  t_ClusteredOptHit->SetBranchAddress("Hit_PE"        , &PDSHit_PE       );
  t_ClusteredOptHit->SetBranchAddress("Hit_OpChannel" , &PDSHit_OpChannel);

  TTree* t_TrueInfo = (TTree*)f_Input->Get("TrueInfo");
  t_TrueInfo->SetBranchAddress("Event",    &EventTrue);
  t_TrueInfo->SetBranchAddress("MarlTime", &MarlTime);
  t_TrueInfo->SetBranchAddress("ENu",      &ENu     );
  t_TrueInfo->SetBranchAddress("ENu_Lep",  &ENu_Lep );
  t_TrueInfo->SetBranchAddress("PosX",     &PosX    );
  t_TrueInfo->SetBranchAddress("PosY",     &PosY    );
  t_TrueInfo->SetBranchAddress("PosZ",     &PosZ    );
  t_TrueInfo->SetBranchAddress("PosT",     &PosT    );
  t_TrueInfo->SetBranchAddress("DirX",     &DirX    );
  t_TrueInfo->SetBranchAddress("DirY",     &DirY    );
  t_TrueInfo->SetBranchAddress("DirZ",     &DirZ    );

  std::map<int,int> map_Event_TrueEntry;
  int nMarleyEvent=0;
  for(int i = 0; i < t_TrueInfo->GetEntries(); i++) {
    t_TrueInfo->GetEntry(i);
    map_Event_TrueEntry[EventTrue] = i;
    nMarleyEvent += (int)ENu->size();
  }
  std::cout << "There were " << nMarleyEvent << " Marley events in " << map_Event_TrueEntry.size() << " Larsoft events." << std::endl;

  std::map<int,std::vector<int> > map_event_entry_opti;

  std::cout << t_ClusteredOptHit->GetEntries() << " optical cluster were saved." << std::endl;

  for(int i=0; i<t_ClusteredOptHit->GetEntries();++i) {
    t_ClusteredOptHit->GetEntry(i);
    if (Config==InputConfig) map_event_entry_opti[Event].push_back(i);
  }

  TCanvas c;
  c.Print("opt_clustering.pdf[");
  TProfile* p_gentype_sign_opti = SetUpTProfileGenType("p_gentype_sign_opti", ";;average number of hits in cluster");
  TProfile* p_gentype_back_opti = SetUpTProfileGenType("p_gentype_back_opti", ";;average number of hits in cluster");
  p_gentype_sign_opti->SetLineColor(kRed);
  p_gentype_back_opti->SetLineColor(kBlue);
  p_gentype_sign_opti->SetLineStyle(1);
  p_gentype_back_opti->SetLineStyle(1);
  p_gentype_sign_opti->SetLineWidth(2);
  p_gentype_back_opti->SetLineWidth(2);


  TH1D* h_ncluster_sign_opti = new TH1D("h_ncluster_sign_opti", ";n clusters;Events", 1500, 0, 1500);
  TH1D* h_ncluster_back_opti = new TH1D("h_ncluster_back_opti", ";n clusters;Events", 1500, 0, 1500);
  h_ncluster_sign_opti->SetLineColor(kRed);
  h_ncluster_back_opti->SetLineColor(kBlue);
  h_ncluster_sign_opti->SetLineStyle(1);
  h_ncluster_back_opti->SetLineStyle(1);
  h_ncluster_sign_opti->SetLineWidth(2);
  h_ncluster_back_opti->SetLineWidth(2);

  TH1D* h_nhit_sign_opti = new TH1D("h_nhit_sign_opti", ";n Hits;Clusters", 100, -0.5, 99.5);
  TH1D* h_nhit_back_opti = new TH1D("h_nhit_back_opti", ";n Hits;Clusters", 100, -0.5, 99.5);
  h_nhit_sign_opti->SetLineColor(kRed);
  h_nhit_back_opti->SetLineColor(kBlue);
  h_nhit_sign_opti->SetLineStyle(1);
  h_nhit_back_opti->SetLineStyle(1);
  h_nhit_sign_opti->SetLineWidth(2);
  h_nhit_back_opti->SetLineWidth(2);

  TEfficiency* eff_enu_opti  = new TEfficiency("eff_enu_opti",  ";E#nu [MeV];Efficiency", 10, 0, 30);

  t_ClusteredOptHit->Project("h_nhit_sign_opti", "NHits", Form("Type==1 && Config==%i",InputConfig));
  t_ClusteredOptHit->Project("h_nhit_back_opti", "NHits", Form("Type==0 && Config==%i",InputConfig));
  std::vector<double> max = {h_nhit_sign_opti->GetMaximum(),
                             h_nhit_back_opti->GetMaximum()};

  h_nhit_sign_opti->SetMaximum((*std::max_element(max.begin(),max.end()))*2);
  h_nhit_sign_opti->SetMinimum(0.1);
  gPad->SetLogy();
  h_nhit_sign_opti->Draw("");
  h_nhit_back_opti->Draw("SAME");
  c.Print("opt_clustering.pdf");

  TH1D* h_npe_sign_opti = new TH1D("h_npe_sign_opti", ";n PEs;Clusters", 20, 0, 20);
  TH1D* h_npe_back_opti = new TH1D("h_npe_back_opti", ";n PEs;Clusters", 20, 0, 20);
  h_npe_sign_opti->SetLineColor(kRed);
  h_npe_back_opti->SetLineColor(kBlue);
  h_npe_sign_opti->SetLineStyle(1);
  h_npe_back_opti->SetLineStyle(1);
  h_npe_sign_opti->SetLineWidth(2);
  h_npe_back_opti->SetLineWidth(2);

  t_ClusteredOptHit->Project("h_npe_sign_opti", "SumPE", Form("Type==1 && Config==%i",InputConfig));
  t_ClusteredOptHit->Project("h_npe_back_opti", "SumPE", Form("Type==0 && Config==%i",InputConfig));
  max = {h_npe_sign_opti->GetMaximum(),
         h_npe_back_opti->GetMaximum()};
  h_npe_sign_opti->SetMaximum((*std::max_element(max.begin(),max.end()))*2);
  h_npe_sign_opti->SetMinimum(0.1);
  gPad->SetLogy();
  h_npe_sign_opti->SetBinContent(h_npe_sign_opti->GetXaxis()->GetNbins(),
                                 h_npe_sign_opti->GetBinContent(h_npe_sign_opti->GetXaxis()->GetNbins()+1));
  h_npe_back_opti->SetBinContent(h_npe_back_opti->GetXaxis()->GetNbins(),
                                 h_npe_back_opti->GetBinContent(h_npe_back_opti->GetXaxis()->GetNbins()+1));
  h_npe_sign_opti->Draw("");
  h_npe_back_opti->Draw("SAME");
  c.Print("opt_clustering.pdf");

  TH1D* h_twidth_sign_opti = new TH1D("h_twidth_sign_opti", ";Time Width [#mus];Clusters", 20,0, 10);
  TH1D* h_twidth_back_opti = new TH1D("h_twidth_back_opti", ";Time Width [#mus];Clusters", 20,0, 10);
  h_twidth_sign_opti->SetLineColor(kRed);
  h_twidth_back_opti->SetLineColor(kBlue);
  h_twidth_sign_opti->SetLineStyle(1);
  h_twidth_back_opti->SetLineStyle(1);
  h_twidth_sign_opti->SetLineWidth(2);
  h_twidth_back_opti->SetLineWidth(2);

  t_ClusteredOptHit->Project("h_twidth_sign_opti", "TimeWidth", Form("Type==1 && Config==%i",InputConfig));
  t_ClusteredOptHit->Project("h_twidth_back_opti", "TimeWidth", Form("Type==0 && Config==%i",InputConfig));
  h_twidth_sign_opti->SetStats(0);
  h_twidth_back_opti->SetStats(0);
  AddOverflow(h_twidth_sign_opti);
  AddOverflow(h_twidth_back_opti);
  max = {h_twidth_sign_opti->GetMaximum(),
         h_twidth_back_opti->GetMaximum()};
  h_twidth_sign_opti->SetMaximum((*std::max_element(max.begin(),max.end()))*2);
  h_twidth_sign_opti->SetMinimum(0.1);
  gPad->SetLogy();
  h_twidth_sign_opti->Draw("");
  h_twidth_back_opti->Draw("SAME");
  c.Print("opt_clustering.pdf");

  TH1D* h_width_sign_opti = new TH1D("h_width_sign_opti", ";Z Width [cm];Clusters", 6, 0., 6.*232.39);
  TH1D* h_width_back_opti = new TH1D("h_width_back_opti", ";Z Width [cm];Clusters", 6, 0., 6.*232.39);
  h_width_sign_opti->SetLineColor(kRed);
  h_width_back_opti->SetLineColor(kBlue);
  h_width_sign_opti->SetLineStyle(1);
  h_width_back_opti->SetLineStyle(1);
  h_width_sign_opti->SetLineWidth(2);
  h_width_back_opti->SetLineWidth(2);

  t_ClusteredOptHit->Project("h_width_sign_opti", "ZWidth", Form("Type==1 && Config==%i",InputConfig));
  t_ClusteredOptHit->Project("h_width_back_opti", "ZWidth", Form("Type==0 && Config==%i",InputConfig));
  max = {h_width_sign_opti->GetMaximum(),
         h_width_back_opti->GetMaximum()};
  h_width_sign_opti->SetMaximum((*std::max_element(max.begin(),max.end()))*2);
  h_width_sign_opti->SetMinimum(0.1);
  gPad->SetLogy();
  h_width_sign_opti->Draw("");
  h_width_back_opti->Draw("SAME");
  c.Print("opt_clustering.pdf");

  TH1D* h_ywidth_sign_opti = new TH1D("h_ywidth_sign_opti", ";Y Width [cm];Clusters", 25, 0., 25.*62.3567);
  TH1D* h_ywidth_back_opti = new TH1D("h_ywidth_back_opti", ";Y Width [cm];Clusters", 25, 0., 25.*62.3567);
  h_ywidth_sign_opti->SetLineColor(kRed);
  h_ywidth_back_opti->SetLineColor(kBlue);
  h_ywidth_sign_opti->SetLineStyle(1);
  h_ywidth_back_opti->SetLineStyle(1);
  h_ywidth_sign_opti->SetLineWidth(2);
  h_ywidth_back_opti->SetLineWidth(2);

  t_ClusteredOptHit->Project("h_ywidth_sign_opti", "YWidth", Form("Type==1 && Config==%i",InputConfig));
  t_ClusteredOptHit->Project("h_ywidth_back_opti", "YWidth", Form("Type==0 && Config==%i",InputConfig));
  max = {h_ywidth_sign_opti->GetMaximum(),
         h_ywidth_back_opti->GetMaximum()};
  h_ywidth_sign_opti->SetMaximum((*std::max_element(max.begin(),max.end()))*2);
  h_ywidth_sign_opti->SetMinimum(0.1);
  h_ywidth_sign_opti->Draw("");
  h_ywidth_back_opti->Draw("SAME");
  c.Print("opt_clustering.pdf");
  int nEventOpti=0;
  int nDetectedSignalEventOpti=0;
  int nBackgroundEventOpti=0;
  std::cout << "map_event_entry_opti " << map_event_entry_opti.size() << std::endl;

  for(auto const& it : map_event_entry_opti) {
    int ncluster_sign=0;
    int ncluster_back=0;
    ++nEventOpti;
    bool SignDetected=false;
    for (auto const& it2 : it.second) {
      std::map<int, int> map_gentype_nhit_sign;
      std::map<int, int> map_gentype_nhit_back;
      t_ClusteredOptHit->GetEntry(it2);
      if (nHitCut!=-1 && PDSHit_GenType->size()<nHitCut) continue;
      if (Type==0) {
        map_gentype_nhit_back = GetMapOfHit(PDSHit_GenType);
        ++ncluster_back;
        ++nBackgroundEventOpti;
      } else {
        SignDetected=true;
        map_gentype_nhit_sign = GetMapOfHit(PDSHit_GenType);
        ++ncluster_sign;
      }
      for (auto const& genhit: map_gentype_nhit_sign) p_gentype_sign_opti->Fill(genhit.first, genhit.second);
      for (auto const& genhit: map_gentype_nhit_back) p_gentype_back_opti->Fill(genhit.first, genhit.second);
    }
    if(SignDetected)
      ++nDetectedSignalEventOpti;
    t_TrueInfo->GetEntry(map_Event_TrueEntry[Event]);
    double Energy_Neutrino = ENu->at(MarleyIndex);
    eff_enu_opti->Fill(SignDetected, Energy_Neutrino*1000.);
    h_ncluster_sign_opti->Fill(ncluster_sign);
    h_ncluster_back_opti->Fill(ncluster_back);
  }
  std::cout << "Here " <<std::endl;

  gPad->SetLogy(false);
  eff_enu_opti->Draw();
  c.Print("opt_clustering.pdf");
  gPad->SetLogy(true);

  max = {p_gentype_sign_opti->GetMaximum(),
         p_gentype_back_opti->GetMaximum()};

  p_gentype_sign_opti->SetMaximum((*std::max_element(max.begin(),max.end()))*2);
  p_gentype_sign_opti->SetMinimum(0.1);
  gPad->SetGridx();
  gPad->SetGridy();
  p_gentype_sign_opti->SetStats(0);
  p_gentype_sign_opti->Draw("E");
  p_gentype_back_opti->Draw("E SAME");
  c.Print("opt_clustering.pdf");
  gPad->SetGridx(false);
  gPad->SetGridy(false);


  max = {h_ncluster_sign_opti->GetMaximum(),
         h_ncluster_back_opti->GetMaximum()};

  h_ncluster_sign_opti->SetMaximum((*std::max_element(max.begin(),max.end()))*2);
  h_ncluster_sign_opti->SetMinimum(0.1);
  gPad->SetLogx();

  auto axis = h_ncluster_sign_opti->GetXaxis();
  axis->SetLimits(0.5,1000);
  h_ncluster_sign_opti->Draw("");
  h_ncluster_back_opti->Draw("SAME");
  c.Print("opt_clustering.pdf");

  gPad->Clear();
  TLegend* leg = new TLegend(0.1, 0.1, 0.9, 0.9);
  leg->AddEntry(p_gentype_sign_opti, "signal optical clusters");
  leg->AddEntry(p_gentype_back_opti, "background optical clusters");
  leg->Draw();
  c.Print("opt_clustering.pdf");
  c.Print("opt_clustering.pdf]");
  return 1;

};
