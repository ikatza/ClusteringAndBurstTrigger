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
  while ( (opt = getopt(argc, argv, "n:")) != -1 ) {  // for each option...
    switch ( opt ) {
    case 'n':
      nHitCut = atoi(optarg);
      break;
    case '?':  // unknown option...
      std::cerr << "Unknown option: '" << char(optopt) << "'!" << std::endl;
      break;
    }
  }

  
  TFile *f_Input = new TFile("/dune/app/users/plasorak/ClusteredDevel.root", "READ");

  int    Cluster       ;
  int    Event         ;
  int    Config        ;
  int    StartChan     ;
  int    EndChan       ;
  int    ChanWidth     ;
  int    NChan         ;
  int    Type          ;
  int    NHits         ;
  float  SumADC        ;
  float  FirstTimeHit  ;
  float  LastTimeHit   ;
  float  TimeWidth     ;
  double ENu           ;
  double E_deposited   ;
  double MC_UnderlyingE;
  double ENu_Lep       ;
  double MarlTime      ;
  double TrClusterPosX ;
  double TrClusterPosY ;
  double TrClusterPosZ ;
  double RecClusterPosX;
  double RecClusterPosY;
  double RecClusterPosZ;
  double pur_Other     ;
  double pur_SNnu      ;
  double pur_APA       ;
  double pur_CPA       ;
  double pur_Ar39      ;
  double pur_Neutron   ;
  double pur_Krypton   ;
  double pur_Polonium  ;
  double pur_Radon     ;
  double pur_Ar42      ;
  double YWidth        ;
  double ZWidth        ;
  double SumPE         ;
  
  std::vector<int>    * HitView  = NULL;
  std::vector<int>    * GenType  = NULL;
  std::vector<int>    * HitChan  = NULL;
  std::vector<double> * HitTime  = NULL;
  std::vector<double> * HitSADC  = NULL;
  std::vector<double> * HitRMS   = NULL;
  std::vector<double> * TrPosX   = NULL; 
  std::vector<double> * TrPosY   = NULL; 
  std::vector<double> * TrPosZ   = NULL;

  TTree* t_Output_triggeredclusteredhits = (TTree*)f_Input->Get("Clustering_triggeredclusteredhits");
  t_Output_triggeredclusteredhits->SetBranchAddress("Cluster",        &Cluster      );
  t_Output_triggeredclusteredhits->SetBranchAddress("Event",          &Event        );
  t_Output_triggeredclusteredhits->SetBranchAddress("Config",         &Config       );
  t_Output_triggeredclusteredhits->SetBranchAddress("StartChan",      &StartChan    );
  t_Output_triggeredclusteredhits->SetBranchAddress("EndChan",        &EndChan      );
  t_Output_triggeredclusteredhits->SetBranchAddress("ChanWidth",      &ChanWidth    );
  t_Output_triggeredclusteredhits->SetBranchAddress("NChan",          &NChan        );
  t_Output_triggeredclusteredhits->SetBranchAddress("Type",           &Type         );
  t_Output_triggeredclusteredhits->SetBranchAddress("NHits",          &NHits        );
  t_Output_triggeredclusteredhits->SetBranchAddress("SumADC",         &SumADC       );
  t_Output_triggeredclusteredhits->SetBranchAddress("FirstTimeHit",   &FirstTimeHit );
  t_Output_triggeredclusteredhits->SetBranchAddress("LastTimeHit",    &LastTimeHit  );
  t_Output_triggeredclusteredhits->SetBranchAddress("TimeWidth",      &TimeWidth    );
  t_Output_triggeredclusteredhits->SetBranchAddress("ENu",            &ENu          );
  t_Output_triggeredclusteredhits->SetBranchAddress("E_deposited",    &E_deposited  );
  t_Output_triggeredclusteredhits->SetBranchAddress("MC_UnderlyingE", &MC_UnderlyingE);
  t_Output_triggeredclusteredhits->SetBranchAddress("ENu_Lep",        &ENu_Lep       );
  t_Output_triggeredclusteredhits->SetBranchAddress("MarlTime",       &MarlTime      );
  t_Output_triggeredclusteredhits->SetBranchAddress("RecClusterPosX", &RecClusterPosX);
  t_Output_triggeredclusteredhits->SetBranchAddress("RecClusterPosY", &RecClusterPosY);
  t_Output_triggeredclusteredhits->SetBranchAddress("RecClusterPosZ", &RecClusterPosZ);
  t_Output_triggeredclusteredhits->SetBranchAddress("TrClusterPosX",  &TrClusterPosX);
  t_Output_triggeredclusteredhits->SetBranchAddress("TrClusterPosY",  &TrClusterPosY);
  t_Output_triggeredclusteredhits->SetBranchAddress("TrClusterPosZ",  &TrClusterPosZ);
  t_Output_triggeredclusteredhits->SetBranchAddress("pur_Other"   ,   &pur_Other   );
  t_Output_triggeredclusteredhits->SetBranchAddress("pur_SNnu"    ,   &pur_SNnu    );
  t_Output_triggeredclusteredhits->SetBranchAddress("pur_APA"     ,   &pur_APA     );
  t_Output_triggeredclusteredhits->SetBranchAddress("pur_CPA"     ,   &pur_CPA     );
  t_Output_triggeredclusteredhits->SetBranchAddress("pur_Ar39"    ,   &pur_Ar39    );
  t_Output_triggeredclusteredhits->SetBranchAddress("pur_Neutron" ,   &pur_Neutron );
  t_Output_triggeredclusteredhits->SetBranchAddress("pur_Krypton" ,   &pur_Krypton );
  t_Output_triggeredclusteredhits->SetBranchAddress("pur_Polonium",   &pur_Polonium);
  t_Output_triggeredclusteredhits->SetBranchAddress("pur_Radon"   ,   &pur_Radon   );
  t_Output_triggeredclusteredhits->SetBranchAddress("pur_Ar42"    ,   &pur_Ar42    );
  t_Output_triggeredclusteredhits->SetBranchAddress("HitView",        &HitView);
  t_Output_triggeredclusteredhits->SetBranchAddress("GenType",        &GenType);
  t_Output_triggeredclusteredhits->SetBranchAddress("HitChan",        &HitChan);
  t_Output_triggeredclusteredhits->SetBranchAddress("HitTime",        &HitTime);
  t_Output_triggeredclusteredhits->SetBranchAddress("HitSADC",        &HitSADC);
  t_Output_triggeredclusteredhits->SetBranchAddress("HitRMS",         &HitRMS );
  t_Output_triggeredclusteredhits->SetBranchAddress("TrPosX",         &TrPosX );
  t_Output_triggeredclusteredhits->SetBranchAddress("TrPosY",         &TrPosY );
  t_Output_triggeredclusteredhits->SetBranchAddress("TrPosZ",         &TrPosZ );
  
  TTree* t_Output_triggeredclusteredopthits = (TTree*)f_Input->Get("Clustering_triggeredclusteredopthits");
  t_Output_triggeredclusteredopthits->SetBranchAddress("Cluster",        &Cluster       );
  t_Output_triggeredclusteredopthits->SetBranchAddress("Config",         &Config        );
  t_Output_triggeredclusteredopthits->SetBranchAddress("Event",          &Event         );
  t_Output_triggeredclusteredopthits->SetBranchAddress("YWidth",         &YWidth        );
  t_Output_triggeredclusteredopthits->SetBranchAddress("ZWidth",         &ZWidth        );
  t_Output_triggeredclusteredopthits->SetBranchAddress("NChan",          &NChan         );
  t_Output_triggeredclusteredopthits->SetBranchAddress("Type",           &Type          );
  t_Output_triggeredclusteredopthits->SetBranchAddress("NHits",          &NHits         );
  t_Output_triggeredclusteredopthits->SetBranchAddress("SumPE",          &SumPE         );
  t_Output_triggeredclusteredopthits->SetBranchAddress("FirstTimeHit",   &FirstTimeHit  );
  t_Output_triggeredclusteredopthits->SetBranchAddress("LastTimeHit",    &LastTimeHit   );
  t_Output_triggeredclusteredopthits->SetBranchAddress("TimeWidth",      &TimeWidth     );
  t_Output_triggeredclusteredopthits->SetBranchAddress("ENu",            &ENu           );
  t_Output_triggeredclusteredopthits->SetBranchAddress("RecClusterPosX", &RecClusterPosX);
  t_Output_triggeredclusteredopthits->SetBranchAddress("RecClusterPosY", &RecClusterPosY);
  t_Output_triggeredclusteredopthits->SetBranchAddress("RecClusterPosZ", &RecClusterPosZ);
  t_Output_triggeredclusteredopthits->SetBranchAddress("GenType",        &GenType);
  t_Output_triggeredclusteredopthits->SetBranchAddress("HitTime",        &HitTime);
  t_Output_triggeredclusteredopthits->SetBranchAddress("HitSADC",        &HitSADC);
  t_Output_triggeredclusteredopthits->SetBranchAddress("HitRMS",         &HitRMS );

  TCanvas c;
  c.Print("opt_clustering.pdf[");
  TProfile* p_gentype_sign_wire = SetUpTProfileGenType("p_gentype_sign_wire", ";;average number of hits in cluster");
  TProfile* p_gentype_back_wire = SetUpTProfileGenType("p_gentype_back_wire", ";;average number of hits in cluster");
  TProfile* p_gentype_sign_opti = SetUpTProfileGenType("p_gentype_sign_opti", ";;average number of hits in cluster");
  TProfile* p_gentype_back_opti = SetUpTProfileGenType("p_gentype_back_opti", ";;average number of hits in cluster");
  p_gentype_sign_wire->SetLineColor(kRed);
  p_gentype_back_wire->SetLineColor(kBlue);
  p_gentype_sign_opti->SetLineColor(kRed);
  p_gentype_back_opti->SetLineColor(kBlue);
  p_gentype_sign_wire->SetLineStyle(1);
  p_gentype_back_wire->SetLineStyle(1);
  p_gentype_sign_opti->SetLineStyle(1);
  p_gentype_back_opti->SetLineStyle(1);
  p_gentype_sign_wire->SetLineWidth(2);
  p_gentype_back_wire->SetLineWidth(2);
  p_gentype_sign_opti->SetLineWidth(2);
  p_gentype_back_opti->SetLineWidth(2);

  std::map<int,std::vector<int> > map_event_entry_wire;
  std::map<int,std::vector<int> > map_event_entry_opti;
  for(int i=0; i<t_Output_triggeredclusteredhits->GetEntries();++i) {
    t_Output_triggeredclusteredhits->GetEntry(i);
    if (Config==5) map_event_entry_wire[Event].push_back(i);
  }

  for(int i=0; i<t_Output_triggeredclusteredopthits->GetEntries();++i) {
    t_Output_triggeredclusteredopthits->GetEntry(i);
    if (Config==5) map_event_entry_opti[Event].push_back(i);
  }

  TH1D* h_ncluster_sign_wire = new TH1D("h_ncluster_sign_wire", ";n clusters;Events", 1500, 0, 1500);
  TH1D* h_ncluster_back_wire = new TH1D("h_ncluster_back_wire", ";n clusters;Events", 1500, 0, 1500);
  TH1D* h_ncluster_sign_opti = new TH1D("h_ncluster_sign_opti", ";n clusters;Events", 1500, 0, 1500);
  TH1D* h_ncluster_back_opti = new TH1D("h_ncluster_back_opti", ";n clusters;Events", 1500, 0, 1500);
  h_ncluster_sign_wire->SetLineColor(kRed);
  h_ncluster_back_wire->SetLineColor(kBlue);
  h_ncluster_sign_opti->SetLineColor(kRed);
  h_ncluster_back_opti->SetLineColor(kBlue);
  h_ncluster_sign_wire->SetLineStyle(1);
  h_ncluster_back_wire->SetLineStyle(1);
  h_ncluster_sign_opti->SetLineStyle(1);
  h_ncluster_back_opti->SetLineStyle(1);
  h_ncluster_sign_wire->SetLineWidth(2);
  h_ncluster_back_wire->SetLineWidth(2);
  h_ncluster_sign_opti->SetLineWidth(2);
  h_ncluster_back_opti->SetLineWidth(2);

  TH1D* h_nhit_sign_wire = new TH1D("h_nhit_sign_wire", ";n Hits;Clusters", 100, -0.5, 99.5);
  TH1D* h_nhit_back_wire = new TH1D("h_nhit_back_wire", ";n Hits;Clusters", 100, -0.5, 99.5);
  TH1D* h_nhit_sign_opti = new TH1D("h_nhit_sign_opti", ";n Hits;Clusters", 100, -0.5, 99.5);
  TH1D* h_nhit_back_opti = new TH1D("h_nhit_back_opti", ";n Hits;Clusters", 100, -0.5, 99.5);
  h_nhit_sign_wire->SetLineColor(kRed);
  h_nhit_back_wire->SetLineColor(kBlue);
  h_nhit_sign_opti->SetLineColor(kRed);
  h_nhit_back_opti->SetLineColor(kBlue);
  h_nhit_sign_wire->SetLineStyle(1);
  h_nhit_back_wire->SetLineStyle(1);
  h_nhit_sign_opti->SetLineStyle(1);
  h_nhit_back_opti->SetLineStyle(1);
  h_nhit_sign_wire->SetLineWidth(2);
  h_nhit_back_wire->SetLineWidth(2);
  h_nhit_sign_opti->SetLineWidth(2);
  h_nhit_back_opti->SetLineWidth(2);

  TEfficiency* eff_enu_opti  = new TEfficiency("eff_enu_opti",  ";E#nu [MeV];Efficiency", 10, 0, 30);
  // TEfficiency* eff_enu_wire  = new TEfficiency("eff_enu_wire",  ";E#nu [MeV];Efficiency");
  // TEfficiency* eff_posx_opti = new TEfficiency("eff_posx_opti", ";Position x [cm];Efficiency");
  // TEfficiency* eff_posy_opti = new TEfficiency("eff_posy_opti", ";Position y [cm];Efficiency");
  // TEfficiency* eff_posz_opti = new TEfficiency("eff_posz_opti", ";Position z [cm];Efficiency");
  // TEfficiency* eff_posx_wire = new TEfficiency("eff_posx_wire", ";Position x [cm];Efficiency");
  // TEfficiency* eff_posy_wire = new TEfficiency("eff_posy_wire", ";Position y [cm];Efficiency");
  // TEfficiency* eff_posz_wire = new TEfficiency("eff_posz_wire", ";Position z [cm];Efficiency");
  // TEfficiency* eff_dirx_opti = new TEfficiency("eff_dirx_opti", ";Direction x;Efficiency");
  // TEfficiency* eff_diry_opti = new TEfficiency("eff_diry_opti", ";Direction y;Efficiency");
  // TEfficiency* eff_dirz_opti = new TEfficiency("eff_dirz_opti", ";Direction z;Efficiency");
  // TEfficiency* eff_dirx_wire = new TEfficiency("eff_dirx_wire", ";Direction x;Efficiency");
  // TEfficiency* eff_diry_wire = new TEfficiency("eff_diry_wire", ";Direction y;Efficiency");
  // TEfficiency* eff_dirz_wire = new TEfficiency("eff_dirz_wire", ";Direction z;Efficiency");


  t_Output_triggeredclusteredopthits->Project("h_nhit_sign_opti", "NHits", "Type==1 && Config==5");
  t_Output_triggeredclusteredopthits->Project("h_nhit_back_opti", "NHits", "Type==0 && Config==5");
  t_Output_triggeredclusteredhits   ->Project("h_nhit_sign_wire", "NHits", "Type==1 && Config==5");
  t_Output_triggeredclusteredhits   ->Project("h_nhit_back_wire", "NHits", "Type==0 && Config==5");
  std::vector<double> max = {h_nhit_sign_wire->GetMaximum(),
                             h_nhit_back_wire->GetMaximum(),
                             h_nhit_sign_opti->GetMaximum(),
                             h_nhit_back_opti->GetMaximum()};
  
  h_nhit_sign_wire->SetMaximum((*std::max_element(max.begin(),max.end()))*2);
  h_nhit_sign_wire->SetMinimum(0.1);
  h_nhit_sign_opti->SetMaximum((*std::max_element(max.begin(),max.end()))*2);
  h_nhit_sign_opti->SetMinimum(0.1);
  gPad->SetLogy();
  // h_nhit_sign_wire->Draw();
  // h_nhit_back_wire->Draw("SAME");
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

  t_Output_triggeredclusteredopthits   ->Project("h_npe_sign_opti", "SumPE", "Type==1 && Config==5");
  t_Output_triggeredclusteredopthits   ->Project("h_npe_back_opti", "SumPE", "Type==0 && Config==5");
  max = {h_npe_sign_opti->GetMaximum(),
         h_npe_back_opti->GetMaximum()};
  h_npe_sign_opti->SetMaximum((*std::max_element(max.begin(),max.end()))*2);
  h_npe_sign_opti->SetMinimum(0.1);
  gPad->SetLogy();
  h_npe_sign_opti->SetBinContent(h_npe_sign_opti->GetXaxis()->GetNbins(),
                                 h_npe_sign_opti->GetBinContent(h_npe_sign_opti->GetXaxis()->GetNbins()+1));
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

  t_Output_triggeredclusteredopthits   ->Project("h_twidth_sign_opti", "TimeWidth", "Type==1 && Config==5");
  t_Output_triggeredclusteredopthits   ->Project("h_twidth_back_opti", "TimeWidth", "Type==0 && Config==5");
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

  t_Output_triggeredclusteredopthits   ->Project("h_width_sign_opti", "ZWidth", "Type==1 && Config==5");
  t_Output_triggeredclusteredopthits   ->Project("h_width_back_opti", "ZWidth", "Type==0 && Config==5");
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

  t_Output_triggeredclusteredopthits->Project("h_ywidth_sign_opti", "YWidth", "Type==1 && Config==5");
  t_Output_triggeredclusteredopthits->Project("h_ywidth_back_opti", "YWidth", "Type==0 && Config==5");
  max = {h_ywidth_sign_opti->GetMaximum(),
         h_ywidth_back_opti->GetMaximum()};
  h_ywidth_sign_opti->SetMaximum((*std::max_element(max.begin(),max.end()))*2);
  h_ywidth_sign_opti->SetMinimum(0.1);
  h_ywidth_sign_opti->Draw("");
  h_ywidth_back_opti->Draw("SAME");
  c.Print("opt_clustering.pdf");
  
  // eff_enu_wire ->Fill();
  // eff_posx_wire->Fill();
  // eff_posy_wire->Fill();
  // eff_posz_wire->Fill();
  // eff_dirx_wire->Fill();
  // eff_diry_wire->Fill();
  // eff_dirz_wire->Fill();

  int nEventWire=0;
  int nDetectedSignalEventWire=0;
  int nBackgroundEventWire=0;
  for(auto const& it : map_event_entry_wire) {
    int ncluster_sign=0;
    int ncluster_back=0;
    ++nEventWire;
    bool SignDetected=false;
    for (auto const& it2 : it.second) {
      std::map<int, int> map_gentype_nhit_sign;
      std::map<int, int> map_gentype_nhit_back;
      t_Output_triggeredclusteredhits->GetEntry(it2);
      if (Type==0) {
        map_gentype_nhit_back = GetMapOfHit(GenType);
        ++ncluster_back;
        ++nBackgroundEventWire;
      } else {
        SignDetected=true;
        map_gentype_nhit_sign = GetMapOfHit(GenType);
        ++ncluster_sign;
      }
      for (auto const& genhit: map_gentype_nhit_sign) p_gentype_sign_wire->Fill(genhit.first, genhit.second);
      for (auto const& genhit: map_gentype_nhit_back) p_gentype_back_wire->Fill(genhit.first, genhit.second);
    }
    if(SignDetected)
      ++nDetectedSignalEventWire;
    h_ncluster_sign_wire->Fill(ncluster_sign);
    h_ncluster_back_wire->Fill(ncluster_back);
  }

  int nEventOpti=0;
  int nDetectedSignalEventOpti=0;
  int nBackgroundEventOpti=0;
  for(auto const& it : map_event_entry_opti) {
    int ncluster_sign=0;
    int ncluster_back=0;
    ++nEventOpti;
    bool SignDetected=false;
    for (auto const& it2 : it.second) {
      std::map<int, int> map_gentype_nhit_sign;
      std::map<int, int> map_gentype_nhit_back;
      t_Output_triggeredclusteredopthits->GetEntry(it2);
      if (nHitCut!=-1 && GenType->size()<nHitCut) continue;
      if (Type==0) {
        map_gentype_nhit_back = GetMapOfHit(GenType);
        ++ncluster_back;
        ++nBackgroundEventOpti;
      } else {
        SignDetected=true;
        map_gentype_nhit_sign = GetMapOfHit(GenType);
        ++ncluster_sign;
      }
      for (auto const& genhit: map_gentype_nhit_sign) p_gentype_sign_opti->Fill(genhit.first, genhit.second);
      for (auto const& genhit: map_gentype_nhit_back) p_gentype_back_opti->Fill(genhit.first, genhit.second);
    }
    if(SignDetected)
      ++nDetectedSignalEventOpti;
    eff_enu_opti->Fill(SignDetected, ENu*1000.);
    // eff_posx_opti->Fill(SignDetected, );
    // eff_posy_opti->Fill(SignDetected, );
    // eff_posz_opti->Fill(SignDetected, );
    // eff_dirx_opti->Fill(SignDetected, );
    // eff_diry_opti->Fill(SignDetected, );
    // eff_dirz_opti->Fill(SignDetected, );
    // std::cout << ncluster_sign << std::endl;
    // std::cout << ncluster_back << std::endl;
    h_ncluster_sign_opti->Fill(ncluster_sign);
    h_ncluster_back_opti->Fill(ncluster_back);
  }

  gPad->SetLogy(false);
  eff_enu_opti->Draw();
  c.Print("opt_clustering.pdf");
  gPad->SetLogy(true);

  max = {//p_gentype_sign_wire->GetMaximum(),
         //p_gentype_back_wire->GetMaximum(),
         p_gentype_sign_opti->GetMaximum(),
         p_gentype_back_opti->GetMaximum()};
  
  p_gentype_sign_wire->SetMaximum((*std::max_element(max.begin(),max.end()))*2);
  p_gentype_sign_wire->SetMinimum(0.1);
  p_gentype_sign_opti->SetMaximum((*std::max_element(max.begin(),max.end()))*2);
  p_gentype_sign_opti->SetMinimum(0.1);
  gPad->SetGridx();
  gPad->SetGridy();
  p_gentype_sign_wire->SetStats(0);
  p_gentype_sign_opti->SetStats(0);
  p_gentype_sign_wire->Draw("E");
  p_gentype_back_wire->Draw("E SAME");
  p_gentype_sign_opti->Draw("E");
  p_gentype_back_opti->Draw("E SAME");
  c.Print("opt_clustering.pdf");
  gPad->SetGridx(false);
  gPad->SetGridy(false);


  max = {// h_ncluster_sign_wire->GetMaximum(),
         // h_ncluster_back_wire->GetMaximum(),
         h_ncluster_sign_opti->GetMaximum(),
         h_ncluster_back_opti->GetMaximum()};

  h_ncluster_sign_wire->SetMaximum((*std::max_element(max.begin(),max.end()))*2);
  h_ncluster_sign_wire->SetMinimum(0.1);
  h_ncluster_sign_opti->SetMaximum((*std::max_element(max.begin(),max.end()))*2);
  h_ncluster_sign_opti->SetMinimum(0.1);
  gPad->SetLogx();

  h_ncluster_sign_wire->Draw("");
  h_ncluster_back_wire->Draw("SAME");
  auto axis = h_ncluster_sign_opti->GetXaxis();
  axis->SetLimits(0.5,1000); 
  h_ncluster_sign_opti->Draw("");
  h_ncluster_back_opti->Draw("SAME");
  c.Print("opt_clustering.pdf");

  std::cout << "nEventWire               "<< nEventWire               << std::endl;
  std::cout << "nDetectedSignalEventWire "<< nDetectedSignalEventWire << std::endl;
  std::cout << "nBackgroundEventWire     "<< nBackgroundEventWire     << std::endl;
  std::cout << "nEventOpti               "<< nEventOpti               << std::endl;
  std::cout << "nDetectedSignalEventOpti "<< nDetectedSignalEventOpti << std::endl;
  std::cout << "nBackgroundEventOpti     "<< nBackgroundEventOpti     << std::endl;
  std::cout << "Wire eff       " << (double)nDetectedSignalEventWire/nEventOpti               << std::endl;
  std::cout << "Wire back rate " << (double)nBackgroundEventWire/nEventOpti/(2.246e-3)/0.12   << " Hz" << std::endl;
  std::cout << "Opti eff       " << (double)nDetectedSignalEventOpti/nEventOpti               << std::endl;
  std::cout << "Opti back rate " << (double)nBackgroundEventOpti/nEventOpti/((3.2e-3)*0.12*3) << " Hz" << std::endl;
  std::ofstream outFile;
  outFile.open(Form("Rates_nhits_%i.txt",nHitCut));
  outFile << nHitCut << "  "
          << (double)nDetectedSignalEventOpti/nEventOpti << "  "
          << (double)nBackgroundEventOpti/nEventOpti/((3.2e-3)*0.12*3)
          << std::endl;
  gPad->Clear();
  TLegend* leg = new TLegend(0.1, 0.1, 0.9, 0.9);
  // leg->AddEntry(p_gentype_sign_wire, "signal wire clusters");
  // leg->AddEntry(p_gentype_back_wire, "background wire clusters");
  leg->AddEntry(p_gentype_sign_opti, "signal optical clusters");
  leg->AddEntry(p_gentype_back_opti, "background optical clusters");
  leg->Draw();
  c.Print("opt_clustering.pdf");
  gPad->SetLogy();
  c.Print("opt_clustering.pdf]");
  return 1;

};
