#include "PhotonTimingDisplay.hh"

PhotonTimingDisplay::PhotonTimingDisplay(std::string F, std::string T): Display(F,T) {
  fPE_Mode = false;
  position = NULL;
  energy   = NULL;
  f_map_gen_th1_timing = GetHistos("timing", ";Time [#mus];n Hits", 140, -4, 10);
};


PhotonTimingDisplay::~PhotonTimingDisplay() {
  ResetAll();
};


void PhotonTimingDisplay::AddArrow(double time, GenType gen) {
  TArrow* ar = new TArrow(time, 0.2, time, 0, 0.05,"|>");
  std::vector<int>vec_colors = getColors(0);
  ar->SetLineWidth(1);
  ar->SetFillColor(vec_colors[gen]);
  ar->SetLineColor(kBlack);
  f_map_gen_arrow_truthtiming[gen].push_back(ar);
};


void PhotonTimingDisplay::ResetAll(){
  for(auto & it : f_map_gen_th1_timing) {
    delete it.second;
    it.second = NULL;
  }
  f_map_gen_th1_timing.clear();
  for(auto & it : f_map_gen_arrow_truthtiming) {
    for(auto & it2 : it.second) {
      delete it2;
      it2 = NULL;
    }
    it.second.clear();
  }
  f_map_gen_arrow_truthtiming.clear();
  delete position;
  delete energy;
  position = NULL;
  energy   = NULL;
};


void PhotonTimingDisplay::PlotAll(){
  c->cd();
  double maxY=-1;
  for (auto const& it : f_map_gen_th1_timing)
    if (maxY < it.second->GetMaximum())
      maxY = it.second->GetMaximum();
    
  for (auto& it : f_map_gen_th1_timing) {
    it.second->SetMaximum(maxY*1.2);
    it.second->SetTitle("Timing");
  }
  for (std::map<GenType,TH1D*>::iterator it = f_map_gen_th1_timing.begin();
       it != f_map_gen_th1_timing.end(); ++it) {
    if (it == f_map_gen_th1_timing.begin()) {
      it->second->Draw("");
    } else {
      it->second->Draw("SAME");
    }
  }
  for (auto & it : f_map_gen_arrow_truthtiming) {
    for (auto& it2 : it.second) {
      it2->SetY1(maxY*0.2);
      it2->Draw();
    }
    
  }

  if (position) {
    position->SetX(0);
    position->SetY(1.2*maxY*0.8);
    position->Draw();
  }
  if (energy) {
    energy->SetX(0);
    energy->SetY(1.2*maxY*0.9);
    energy->Draw();
  }
    
};


void PhotonTimingDisplay::DisplayEvent(int event, int type, double time) {

  int requestedEvent = event;
  if (requestedEvent<0) requestedEvent=0;
  c->Print("PhotonTimingDisplay.pdf[");
  bool  found=false;

  int tries=0;
  while(!found) {
    
    if (position) {
      delete position;
      position=NULL;
    }
    if (energy) {
      delete energy;
      energy=NULL;
    }
  
    int nEvent = std::min(requestedEvent, (int)im.GetEntries());
    im.GetEntry(nEvent);
    std::cout << "Plotting event number " << nEvent << std::endl;
    if (type == -1) {
      found=true;
      for(size_t it=0; it<im.PDS_OpHit_OpChannel->size(); ++it) {
        double weight=1;
        if (fPE_Mode) weight=(*im.PDS_OpHit_PE)[it];
        f_map_gen_th1_timing[ConvertIntToGenType((*im.PDS_OpHit_True_GenType)[it])]->Fill((*im.PDS_OpHit_PeakTime)[it]-time, weight);
      }
      position = new TLatex(0,0,Form("X: %.1f   Y: %.1f   Z: %.1f",     (*im.True_VertX)[0],(*im.True_VertY)[0],(*im.True_VertZ)[0]));
      energy   = new TLatex(0,0,Form("#nu energy: %.1f MeV  nHit: %.0f",(*im.True_ENu  )[0]*1000., (double)f_map_gen_th1_timing[kSNnu]->GetEntries()));
      for(size_t it=0; it<im.True_Bck_Time->size(); ++it) {
        AddArrow((*im.True_Bck_Time)[it]/1000.-time, ConvertIntToGenType((*im.True_Bck_Mode)[it]));
      }
    } else {
      GenType gen = ConvertIntToGenType(type);

      double timeoffset=0;
      for(size_t it=0; it<im.PDS_OpHit_OpChannel->size(); ++it) {
        if ((*im.PDS_OpHit_True_GenType)[it] == type) {
          found=true;
          timeoffset = (*im.PDS_OpHit_PeakTime)[it];
          if(gen==kSNnu) { energy = new TLatex(0,0,Form("energy: %.1f MeV",(*im.True_ENu)[0]*1000.)); }
          else           { energy = new TLatex(0,0,Form("energy: %.1f MeV",(*im.PDS_OpHit_True_Energy)[0]*1000.)); }
        }
      }
      if (!found) {
        std::cerr << "Didn't find a event of requested truth type in this event, continuing..." << std::endl;
        timeoffset = 0;
      } else {
        for(size_t it=0; it<im.PDS_OpHit_OpChannel->size(); ++it) {
          double weight=1;
          if (fPE_Mode) weight=(*im.PDS_OpHit_PE)[it];
          f_map_gen_th1_timing[ConvertIntToGenType((*im.PDS_OpHit_True_GenType)[it])]->Fill((*im.PDS_OpHit_PeakTime)[it]-time-timeoffset, weight);
        }
        if (found)
          energy->SetText(0,0,(std::string(energy->GetTitle())+Form("  nHit: %.0f",(double)f_map_gen_th1_timing[gen]->GetEntries())).c_str());
        if(gen == kSNnu)
          position = new TLatex(0,0,Form("X: %.1f   Y: %.1f   Z: %.1f",     (*im.True_VertX)[0],(*im.True_VertY)[0],(*im.True_VertZ)[0]));
        for(size_t it=0; it<im.True_Bck_Time->size(); ++it) {
          AddArrow((*im.True_Bck_Time)[it]/1000.-time-timeoffset, ConvertIntToGenType((*im.True_Bck_Mode)[it]));
        }
      }
    }
    
    if (found) {
      PlotAll();
      c->Print("PhotonTimingDisplay.pdf");
    }
    tries++;
    if(tries>(int)im.GetEntries())
      break;
    //ResetAll();
    requestedEvent=nEvent+1;
  }
  c->Print("PhotonTimingDisplay.pdf]");
  
};
