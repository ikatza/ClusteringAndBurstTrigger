#include "EventInspector.hh"
#include "TMath.h"



class NeutronInspector: public EventInspector{

public:
  TH1D h_timing_relat    ;
  TH1D h_timing_start    ;  
  TH1D h_timing_end      ;  
  TH1D h_distance_relat  ;
  TH1D h_distance_start_x;
  TH1D h_distance_start_y;
  TH1D h_distance_start_z;
  TH1D h_distance_end_x  ;
  TH1D h_distance_end_y  ;
  TH1D h_distance_end_z  ;
  TH1D h_energy          ;
  TH2D h2_energy         ;

  NeutronInspector(int argc, char** argv): EventInspector(argc, argv) {
    SetOutFileName("Neutron_Info.pdf");
    int nbins=100;
    double bins[101];
    for(int i=0; i<101; i++){
      bins[i]=TMath::Power(1.17,i);
    }
    double bins2[101];
    for(int i=0; i<101; i++){
      bins2[100-i]=2.*TMath::Power(0.87,i);
    }
    double bins3[101];
    for(int i=0; i<101; i++){
      bins3[100-i]=10*TMath::Power(0.87,i);
    }
    h_timing_relat     = TH1D("h_timing_relat",     ";Time Diff [#mus];Event",   nbins,   bins2);
    h_timing_start     = TH1D("h_timing_start",     ";Time Start [#mus];Event",  nbins,   bins);  
    h_timing_end       = TH1D("h_timing_end",       ";Time End [#mus];Event",    nbins,   bins);  
    h_distance_relat   = TH1D("h_distance_relat",   ";Distance Diff [cm];Event", 100,     bins3);
    h_distance_start_x = TH1D("h_distance_start_x", ";start x [cm];Event",       100, -1000, 1000);
    h_distance_start_y = TH1D("h_distance_start_y", ";start y [cm];Event",       100, -2000, 2000);
    h_distance_start_z = TH1D("h_distance_start_z", ";start z [cm];Event",       100,     0, 2000);
    h_distance_end_x   = TH1D("h_distance_end_x",   ";end x [cm];Event",         100, -1000, 1000);
    h_distance_end_y   = TH1D("h_distance_end_y",   ";end y [cm];Event",         100, -2000, 2000);
    h_distance_end_z   = TH1D("h_distance_end_z",   ";end z [cm];Event",         100,     0, 2000);
    h_energy   = TH1D("energy", ";Neutron KE [MeV];Event", 20, 0, 5);
    h2_energy  = TH2D("energy2d", ";E#nu [MeV];Neutron KE [MeV]", 25, 10, 35, 20, 0, 5);
    h_energy  .SetLineColor(kBlack);
    h_energy  .SetLineWidth(2);

    h_timing_relat    .SetLineColor(kBlack);
    h_timing_start    .SetLineColor(kBlack);
    h_timing_end      .SetLineColor(kBlack);
    h_distance_relat  .SetLineColor(kBlack);
    h_distance_start_x.SetLineColor(kBlack);
    h_distance_start_y.SetLineColor(kBlack);
    h_distance_start_z.SetLineColor(kBlack);
    h_distance_end_x  .SetLineColor(kBlack);
    h_distance_end_y  .SetLineColor(kBlack);
    h_distance_end_z  .SetLineColor(kBlack);
    h_timing_relat    .SetLineWidth(2);
    h_timing_start    .SetLineWidth(2);
    h_timing_end      .SetLineWidth(2);
    h_distance_relat  .SetLineWidth(2);
    h_distance_start_x.SetLineWidth(2);
    h_distance_start_y.SetLineWidth(2);
    h_distance_start_z.SetLineWidth(2);
    h_distance_end_x  .SetLineWidth(2);
    h_distance_end_y  .SetLineWidth(2);
    h_distance_end_z  .SetLineWidth(2);
  };
  
  ~NeutronInspector(){
  };

  
  void SaveHistos(){
    AddOverflow(h_timing_relat    );
    AddOverflow(h_timing_start    );
    AddOverflow(h_timing_end      );
    AddOverflow(h_distance_relat  );
    AddOverflow(h_distance_start_x);
    AddOverflow(h_distance_start_y);
    AddOverflow(h_distance_start_z);
    AddOverflow(h_distance_end_x  );
    AddOverflow(h_distance_end_y  );
    AddOverflow(h_distance_end_z  );

    gPad->SetLogx();
    h_timing_relat    .Draw(); AddPlot();
    h_timing_start    .Draw(); AddPlot();
    h_timing_end      .Draw(); AddPlot();
    h_distance_relat  .Draw(); AddPlot();
    gPad->SetLogx(false);
    h_distance_start_x.Draw(); AddPlot();
    h_distance_start_y.Draw(); AddPlot();
    h_distance_start_z.Draw(); AddPlot();
    h_distance_end_x  .Draw(); AddPlot();
    h_distance_end_y  .Draw(); AddPlot();
    h_distance_end_z  .Draw(); AddPlot(); 

    h_energy.Draw();
    AddPlot();
    h2_energy.Draw("COLZ");
    AddPlot();
  };

private:
  bool Selected() { return true; };
  void ExecuteInLoop();
  double GetDistance(double x, double y, double z) { return sqrt(x*x + y*y + z*z);};
  
};

void NeutronInspector::ExecuteInLoop(){
  double ENu     = im.True_ENu    ->back()*1000.;

  for(size_t part=0; part<im.True_Bck_Mode->size(); ++part) {
    GenType gen = ConvertIntToGenType((*im.True_Bck_Mode)[part]);
    if (gen == kSNnu) {
    //if (true) {
      if ((*im.True_Bck_PDG)[part] == 2112) {
        //nNeutron++;
        double d = 0;
        // GetDistance((*im.True_Bck_VertX)[part]-(*im.True_Bck_EndVertX)[part],
        //             (*im.True_Bck_VertY)[part]-(*im.True_Bck_EndVertY)[part],
        //             (*im.True_Bck_VertZ)[part]-(*im.True_Bck_EndVertZ)[part]);
        //h_timing_relat    .Fill(-(*im.True_Bck_Time)[part]+(*im.True_Bck_EndTime)[part]);
        h_timing_start    .Fill((*im.True_Bck_Time)[part]);
        //h_timing_end      .Fill((*im.True_Bck_EndTime)[part]);
        h_distance_relat  .Fill(d);
        h_distance_start_x.Fill((*im.True_Bck_VertX)[part]);
        h_distance_start_y.Fill((*im.True_Bck_VertY)[part]);
        h_distance_start_z.Fill((*im.True_Bck_VertZ)[part]);
        // h_distance_end_x  .Fill((*im.True_Bck_EndVertX)[part]);
        // h_distance_end_y  .Fill((*im.True_Bck_EndVertY)[part]);
        // h_distance_end_z  .Fill((*im.True_Bck_EndVertZ)[part]);
        //h_timing  .Fill(VertexT);//abs((*im.True_Bck_TimeAll)[part]-(*im.True_Bck_EndTimeAll)[part]));
        // h_distance.Fill(GetDistance(VertX - (*im.True_Bck_VertXAll)[part],//-(*im.True_Bck_EndVertXAll)[part],
        //                             VertY - (*im.True_Bck_VertYAll)[part],//-(*im.True_Bck_EndVertYAll)[part],
        //                             VertZ - (*im.True_Bck_VertZAll)[part]//-(*im.True_Bck_EndVertZAll)[part]
        //                   ));
       
        h_energy  .Fill((*im.True_Bck_EnergyAll)[part]*1000.);
        h2_energy .Fill(ENu,(*im.True_Bck_EnergyAll)[part]*1000.);
      }
    }
  }
};

int main(int argc, char** argv) {

  NeutronInspector in(argc, argv);
  in.Loop();
  in.SaveHistos();
  return 1;

};
