#include "PhotonWirePositionDisplay.hh"

void PhotonWirePositionDisplay::CreateGraphHit(std::map<int,TGraph*>&map_gr_wire, int hit) {
  int chan = (*im.Hit_Chan)[hit];
  if (!map_gr_wire[chan]) {
    map_gr_wire[chan] = new TGraph(2);
    map_gr_wire[chan]->SetPoint(0, (*im.Hit_Z_start)[hit], (*im.Hit_Y_start)[hit]);
    map_gr_wire[chan]->SetPoint(1, (*im.Hit_Z_end  )[hit], (*im.Hit_Y_end  )[hit]);
    map_gr_wire[chan]->SetLineColor(kViolet);
    map_gr_wire[chan]->SetLineWidth(1);
  } else {
    map_gr_wire[chan]->SetLineWidth(1+map_gr_wire[chan]->GetLineWidth());
  }

}

std::map<GenType, TGraph*> PhotonWirePositionDisplay::SetUpGraph(const std::string Title) {

  std::map<GenType, TGraph*> map_gr;
  std::vector<int>vec_colors = getColors(0);
  for (auto const&  it: h.GenName) {
    map_gr[it.first] = new TGraph(0);
    map_gr[it.first]->SetTitle((Title+it.second).c_str());
    map_gr[it.first]->SetName ((Title+h.ShortGenName[it.first]).c_str());
    map_gr[it.first]->SetMarkerStyle(kStar);
    map_gr[it.first]->SetMarkerSize(3);
    map_gr[it.first]->SetMarkerColor(kRed);

  }
  return map_gr;

};


std::map<GenType, TH2D*> PhotonWirePositionDisplay::SetUpHistos2D(const std::string Title,
                                                                  const std::string XAxis, const std::string YAxis,
                                                                  const int nbinX, const double minX, const double maxX,
                                                                  const int nbinY, const double minY, const double maxY) {
  std::map<GenType, TH2D*> map_h2;
  std::vector<int>vec_colors = getColors(0);
  for (auto const&  it: h.GenName) {
    map_h2[it.first] = new TH2D(("gr_"+it.second).c_str(), (it.second+";"+XAxis+";"+YAxis).c_str(),
                                nbinX, minX, maxX, nbinY, minY, maxY);
    map_h2[it.first]->SetStats(0);
  }
  return map_h2;
};




PhotonWirePositionDisplay::PhotonWirePositionDisplay(std::string F, std::string T):
  Display(F,T),
  NHit           (NULL),
  NHitTot        (NULL),
  NeutrinoEnergy (NULL),
  VertexXPosition(NULL) {
  map_h2_BackTracked = SetUpGraph("Position");
  map_h2_PosHit      = Get2DHistos("Diplay",";Z Position [cm];Y Position [cm]",
                                   6, 0, 1392,
                                   100, -800, 800);
  
}



void PhotonWirePositionDisplay::DisplayEvent(int event, int type){
  im.GetEntry(event);

  for (size_t hit=0; hit<im.PDS_OpHit_OpChannel->size(); ++hit) {
    GenType gen = ConvertIntToGenType((*im.PDS_OpHit_True_GenType)[hit]);
    map_h2_PosHit  [gen]->Fill((*im.PDS_OpHit_Z        )[hit], (*im.PDS_OpHit_Y)[hit], (*im.PDS_OpHit_PE)[hit]);
    map_h2_PosHit  [kAll]->Fill((*im.PDS_OpHit_Z        )[hit], (*im.PDS_OpHit_Y)[hit], (*im.PDS_OpHit_PE)[hit]);
  }

  for (size_t hit=0; hit<im.Hit_Chan->size(); ++hit) {
    GenType gen =  ConvertIntToGenType((*im.Hit_True_GenType)[hit]);
    CreateGraphHit(map_map_gr_wire[gen],  hit);
    CreateGraphHit(map_map_gr_wire[kAll], hit);
  }

  for (size_t iPart=0; iPart<im.True_Bck_VertY->size(); ++iPart) {
    GenType gen = ConvertIntToGenType((*im.True_Bck_Mode)[iPart]);
    map_h2_BackTracked[gen] ->Set(map_h2_BackTracked[gen]->GetN()+1);
    map_h2_BackTracked[gen] ->SetPoint(map_h2_BackTracked[gen]->GetN()-1,(*im.True_Bck_VertZ)[iPart],(*im.True_Bck_VertY)[iPart]);
    map_h2_BackTracked[kAll]->Set(map_h2_BackTracked[kAll]->GetN()+1);
    map_h2_BackTracked[kAll]->SetPoint(map_h2_BackTracked[kAll]->GetN()-1,(*im.True_Bck_VertZ)[iPart],(*im.True_Bck_VertY)[iPart]);
  }
  NHit            = new TLatex(50,  700,Form("N opt hits: %.0f    N wire hits: %d",
                                             map_h2_PosHit[kSNnu]->GetEntries(),
                                             (int)map_map_gr_wire[kSNnu].size()));
  NHitTot         = new TLatex(50,  700,Form("N opt hits: %.0f    N wire hits: %d",
                                             map_h2_PosHit[kAll]->GetEntries(),
                                             (int)map_map_gr_wire[kAll].size()));
  NeutrinoEnergy  = new TLatex(50, -700,Form("#nu energy: %.1f MeV",(*im.True_ENu  )[0]*1000.));
  VertexXPosition = new TLatex(800,-700,Form("X Position: %.0f cm" ,(*im.True_VertX)[0]));
  c->Print("OpHitTrigger_results_evdisplay.pdf[");
  for (auto const& it : map_h2_PosHit) {
    it.second->Draw("COLZ");

    for (auto const& it_gr : map_map_gr_wire[it.first])
      it_gr.second->Draw("LP");

    if (it.first != kAr39 && it.first != kAll && it.first != kKrypton)
      map_h2_BackTracked[it.first]->Draw("P");
    if(it.first == kSNnu || it.first == kAll) {
      NeutrinoEnergy ->Draw();
      VertexXPosition->Draw();
    }
    if (it.first == kSNnu){
      NHit   ->Draw();
    }
    if (it.first == kAll ) NHitTot->Draw();
        
    c->Print("OpHitTrigger_results_evdisplay.pdf");
  }
  c->Print("OpHitTrigger_results_evdisplay.pdf]");

}

PhotonWirePositionDisplay::~PhotonWirePositionDisplay(){

  for (auto it: map_map_gr_wire) {
    for (auto it2: it.second) {
      delete it2.second;
      it2.second = NULL;
    }
    it.second.clear();
  }
  map_map_gr_wire.clear();
  
  std::map<GenType,TGraph*> map_h2_BackTracked;
  for (auto it: map_h2_BackTracked) {
    delete it.second;
    it.second = NULL;
  }
  map_h2_BackTracked.clear();
  
  std::map<GenType,TH2D*> map_h2_PosHit;
  for (auto it: map_h2_PosHit) {
    delete it.second;
    it.second=NULL;
  }
  map_h2_PosHit.clear();
  
  delete NHit           ; NHit            = NULL;
  delete NHitTot        ; NHitTot         = NULL;
  delete NeutrinoEnergy ; NeutrinoEnergy  = NULL;
  delete VertexXPosition; VertexXPosition = NULL;
}
