#include "EventInspector.hh"



class PDSTimingInspector: public EventInspector{

public:
  
  TH1D h_timing;
  TH1D h_timing_relat;
  TH1D h_timing_truth;
  TH1D h_timing_truth_electron;
  TH1D h_timing_truth_photon;
  TH1D h_timing_integral;
  TH1D h_timing_truth_integral;
  TH1D h_timing_energy_weighted_truth;
  TH1D h_timing_energy_weighted_truth_photon;
  TH1D h_timing_energy_weighted_truth_electron;
  TH1D h_timing_energy_weighted_integral;
  TH1D h_energy_truth;
  TH1D h_energy_truth_electron;
  TH1D h_energy_truth_photon;
  TH1D h_PE;
  TH1D h_hit;

  
  PDSTimingInspector(int argc, char** argv): EventInspector(argc, argv) {
    SetOutFileName("PDS_Timing.pdf");
    h_timing                                = TH1D("h_timing",                                ";Time [#mus];N optical hits or N particles [a.u.]",              100, 0, 10);
    h_timing_relat                          = TH1D("h_timing_relat",                          ";Time [#mus];N optical hits",                                    100, 0,  5);
    h_timing_truth                          = TH1D("h_timing_truth",                          ";Time [#mus];N particles [a.u.]",                                100, 0, 10);
    h_timing_truth_electron                 = TH1D("h_timing_truth_electron",                 ";Time [#mus];N particles [a.u.]",                                100, 0, 10);
    h_timing_truth_photon                   = TH1D("h_timing_truth_photon",                   ";Time [#mus];N particles [a.u.]",                                100, 0, 10);
    h_timing_integral                       = TH1D("h_timing_integral",                       ";Time [#mus];Fraction of optical hits or particle or of energy", 100, 0, 10);
    h_timing_truth_integral                 = TH1D("h_timing_truth_integral",                 ";Time [#mus];N particles [a.u.]",                                100, 0, 10);
    h_timing_energy_weighted_truth          = TH1D("h_timing_energy_weighted_truth",          ";Time [#mus];Energy [a.u.]",                                     100, 0, 10);
    h_timing_energy_weighted_truth_electron = TH1D("h_timing_energy_weighted_truth_electron", ";Time [#mus];Energy [a.u.]",                                     100, 0, 10);
    h_timing_energy_weighted_truth_photon   = TH1D("h_timing_energy_weighted_truth_photon",   ";Time [#mus];Energy [a.u.]",                                     100, 0, 10);
    h_timing_energy_weighted_integral       = TH1D("h_timing_energy_weighted_integral",       ";Time [#mus];Fraction of energy [a.u]",                          100, 0, 10);
    h_energy_truth                          = TH1D("h_energy_truth",                          ";True energy [MeV];N optical hits [a.u]", 100, 0, 3);
    h_energy_truth_electron                 = TH1D("h_energy_truth_electron",                 ";True energy [MeV];N optical hits [a.u]", 100, 0, 3);
    h_energy_truth_photon                   = TH1D("h_energy_truth_photon",                   ";True energy [MeV];N optical hits [a.u]", 100, 0, 3);
    h_PE                                    = TH1D("h_pe",     ";PE;N optical hits", 100, 0, 10);
    h_hit                                   = TH1D("h_hit",    ";nhit;N optical hits", 100, -0.5, 99.5);

    h_timing                               .Sumw2();
    h_timing_relat                         .Sumw2();
    h_timing_truth                         .Sumw2();
    h_timing_truth_electron                .Sumw2();
    h_timing_truth_photon                  .Sumw2();
    h_timing_integral                      .Sumw2();
    h_timing_truth_integral                .Sumw2();
    h_timing_energy_weighted_truth         .Sumw2();
    h_timing_energy_weighted_truth_electron.Sumw2();
    h_timing_energy_weighted_truth_photon  .Sumw2();
    h_timing_energy_weighted_integral      .Sumw2();
    h_energy_truth                         .Sumw2();
    h_energy_truth_electron                .Sumw2();
    h_energy_truth_photon                  .Sumw2();
    h_PE                                   .Sumw2();
    h_hit                                  .Sumw2();


  };
  
  ~PDSTimingInspector(){};

private:
  bool Selected() {
    return ((*im.True_VertZ)[0] > 920 ||(*im.True_VertZ)[0] > 460);
  };
  
  void ExecuteInLoop() {
    std::vector<std::pair<int,int> > trackIDs;
    int nhit=0;
    for (size_t truePart=0; truePart<im.True_Bck_IDAll->size(); ++truePart) {
      GenType gen = ConvertIntToGenType((*im.True_Bck_ModeAll)[truePart]);
      if (gen == kSNnu) {
        double energy=(*im.True_Bck_EnergyAll)[truePart]*1000.;
        if(energy<0) {
          std::cout << "Negative energy ?!? " << energy <<std::endl;;
          std::cout << "Time " << (*im.True_Bck_TimeAll)[truePart] <<std::endl;;
          std::cout << "PDG  " << (*im.True_Bck_PDGAll)[truePart] <<std::endl;;
          energy=0;
        }
        if((*im.True_Bck_PDGAll)[truePart]>5000) {
          energy=0;
        }
        h_timing_truth                .Fill((*im.True_Bck_TimeAll)[truePart]);
        h_timing_energy_weighted_truth.Fill((*im.True_Bck_TimeAll)[truePart], energy);
        h_energy_truth                .Fill(energy);
        if ((*im.True_Bck_PDGAll)[truePart] == 11) {
          h_timing_truth_electron                .Fill((*im.True_Bck_TimeAll)[truePart]);
          h_timing_energy_weighted_truth_electron.Fill((*im.True_Bck_TimeAll)[truePart], energy);
          h_energy_truth_electron                .Fill(energy);
        }
        if ((*im.True_Bck_PDGAll)[truePart] == 22) {
          h_timing_truth_photon                .Fill((*im.True_Bck_TimeAll)[truePart]);
          h_timing_energy_weighted_truth_photon.Fill((*im.True_Bck_TimeAll)[truePart], energy);
          h_energy_truth_photon                .Fill(energy);
        }
        trackIDs.push_back(std::make_pair((*im.True_Bck_IDAll)[truePart],(*im.True_Bck_PDGAll)[truePart]));
      }
    }
    
    if (!trackIDs.empty()) {
      std::vector<double> times;
      std::vector<int> indices;
      for (size_t i=0; i<trackIDs.size() && times.size()<2; ++i) {
        int TrackID  = trackIDs[i].first;
       
        for (size_t index=0; index<im.PDS_OpHit_True_TrackID->size(); ++index) {
          if ((*im.PDS_OpHit_True_TrackID)[index] == TrackID) {
            h_PE.Fill((*im.PDS_OpHit_PE)[index]);
            nhit++;
            indices.push_back(index);
          }
        }
      }
      h_hit.Fill(nhit);

      for (auto it : indices)
        times.push_back((*im.PDS_OpHit_PeakTimeAbs)[it]);
      
      if (times.size()>2) {
        std::sort(times.begin(), times.end());
        for (auto time : times) {
          //time -= (*times.begin());
          h_timing.Fill(time);
          h_timing_relat.Fill(time - (*times.begin()));
        }
      }
    }
  };
public:
  void PlotHistos(TH1D& h_1, TH1D& h_2,TH1D& h_3, TH1D& h_4){
    AddOverflow(h_1);
    AddOverflow(h_2);
    AddOverflow(h_3);
    AddOverflow(h_4);
    h_1.SetLineWidth(2);
    h_2.SetLineWidth(2);
    h_3.SetLineWidth(2);
    h_4.SetLineWidth(2);
    h_3.Scale(1/h_2.Integral());
    h_4.Scale(1/h_2.Integral());
    h_2.Scale(1/h_2.Integral());
    h_1.Scale(1/h_1.Integral());
    h_2.SetLineColor(kRed);
    h_3.SetLineColor(kGreen);
    h_4.SetLineColor(kYellow);
    h_1.SetStats(0);
    std::vector<double> max = {h_1.GetMaximum(),
                               h_2.GetMaximum(),
                               h_3.GetMaximum(),
                               h_4.GetMaximum()};
    std::vector<double> min = {h_1.GetMinimum(),
                               h_2.GetMinimum(),
                               h_3.GetMinimum(),
                               h_4.GetMinimum()};
    if(gPad->GetLogy()>0) {
      for(auto& it:min){
        if(it==0) it=100;
      }
      h_1.SetMinimum(0.00001);
      h_2.SetMinimum(0.00001);
      h_3.SetMinimum(0.00001);
      h_4.SetMinimum(0.00001);
    } else {
      h_1.SetMinimum(0);
      h_2.SetMinimum(0);
      h_3.SetMinimum(0);
      h_4.SetMinimum(0);
    }
    h_1.SetMaximum((*std::max_element(max.begin(),max.end()))*1.2);
    h_2.SetMaximum((*std::max_element(max.begin(),max.end()))*1.2);
    h_3.SetMaximum((*std::max_element(max.begin(),max.end()))*1.2);
    h_4.SetMaximum((*std::max_element(max.begin(),max.end()))*1.2);
    h_1.Draw("HIST");
    h_2.Draw("HIST SAME");
    h_3.Draw("HIST SAME");
    h_4.Draw("HIST SAME");
  }

};

int main(int argc, char** argv) {

  PDSTimingInspector in(argc, argv);
  in.Loop();

  AddOverflow(in.h_timing_relat);
  in.h_timing_relat.SetLineWidth(2);
  in.h_timing_relat.Scale(1/in.h_timing_relat.Integral());
  in.h_timing_relat.Draw("HIST");
  in.AddPlot();
  
  in.h_timing.GetYaxis()->SetTitle("N optical hits or energy [a.u.]");
  gPad->SetLogy(true);

  in.PlotHistos(in.h_timing,
                in.h_timing_energy_weighted_truth,
                in.h_timing_energy_weighted_truth_electron,
                in.h_timing_energy_weighted_truth_photon);
  in.AddPlot();
  
  gPad->SetLogy(false);
  double integral=in.h_timing.Integral();
  for(int i=0; i<=in.h_timing.GetXaxis()->GetNbins();++i)
    in.h_timing_integral.SetBinContent(i, in.h_timing.Integral(0,i)/integral);

  integral=in.h_timing_truth.Integral();
  for(int i=0; i<=in.h_timing_truth_integral.GetXaxis()->GetNbins();++i) {
    in.h_timing_truth_integral.SetBinContent(i, in.h_timing_truth.Integral(0,i)/integral);
  }
  
  integral=in.h_timing_energy_weighted_truth.Integral();
  for(int i=0; i<=in.h_timing_energy_weighted_integral.GetXaxis()->GetNbins();++i)
    in.h_timing_energy_weighted_integral.SetBinContent(i, in.h_timing_energy_weighted_truth.Integral(0,i)/integral);

  in.h_timing_integral                .SetLineWidth(2);
  in.h_timing_truth_integral          .SetLineWidth(2);
  in.h_timing_energy_weighted_integral.SetLineWidth(2);

  in.h_timing_integral                .SetLineColor(kBlue);
  in.h_timing_truth_integral          .SetLineColor(kRed);  
  in.h_timing_energy_weighted_integral.SetLineColor(kRed);

  in.h_timing_energy_weighted_integral.SetLineStyle(2);

  in.h_timing_integral.SetStats(0);
  in.h_timing_integral.SetMinimum(0);

  in.h_timing_integral                .Draw("HIST");
  in.h_timing_truth_integral          .Draw("HIST SAME");
  in.h_timing_energy_weighted_integral.Draw("HIST SAME");
  in.AddPlot();
  
  in.h_timing_integral.GetYaxis()->SetTitle("Fraction of hits");
  in.h_timing_integral                .Draw("HIST");
  in.AddPlot();
    
  AddOverflow(in.h_energy_truth         );
  AddOverflow(in.h_energy_truth_electron);
  AddOverflow(in.h_energy_truth_photon  );
  in.h_energy_truth         .SetLineWidth(2);
  in.h_energy_truth_electron.SetLineWidth(2);
  in.h_energy_truth_photon  .SetLineWidth(2);

  in.h_energy_truth_electron.Scale(1/in.h_energy_truth.Integral());
  in.h_energy_truth_photon  .Scale(1/in.h_energy_truth.Integral());
  in.h_energy_truth         .Scale(1/in.h_energy_truth.Integral());

  in.h_energy_truth         .SetLineColor(kRed);
  in.h_energy_truth_electron.SetLineColor(kGreen);
  in.h_energy_truth_photon  .SetLineColor(kYellow);

  in.h_energy_truth.SetStats(0);
  in.h_energy_truth.SetMinimum(0);
  std::vector<double> max = {in.h_energy_truth         .GetMaximum(),
                             in.h_energy_truth_electron.GetMaximum(),
                             in.h_energy_truth_photon  .GetMaximum()};

  in.h_energy_truth.SetMaximum((*std::max_element(max.begin(),max.end()))*1.2);
  in.h_energy_truth         .Draw("HIST");
  in.h_energy_truth_electron.Draw("HIST SAME");
  in.h_energy_truth_photon  .Draw("HIST SAME");
  in.AddPlot();

  gPad->SetLogx(false);
  in.h_PE .Draw();
  in.AddPlot();
  in.h_hit.Draw();
  in.AddPlot();
   
  return 0;
}
