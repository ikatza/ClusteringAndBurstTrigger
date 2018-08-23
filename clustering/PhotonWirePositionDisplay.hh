#ifndef PHOTONWIREPOSITIONDISPLAY_HH
#define PHOTONWIREPOSITIONDISPLAY_HH

#include "Display.hh"
#include "TArrow.h"
#include "TLatex.h"


class PhotonWirePositionDisplay: public Display {

public:
  PhotonWirePositionDisplay(std::string, std::string);
  void DisplayEvent(const int event=-1, const int type=-1);
  ~PhotonWirePositionDisplay();
  void CreateGraphHit(std::map<int,TGraph*>&, int);
  std::map<GenType, TGraph*> SetUpGraph(const std::string);
  std::map<GenType, TH2D*> SetUpHistos2D(const std::string,
                                         const std::string, const std::string,
                                         const int, const double, const double,
                                         const int, const double, const double);

private:
  std::map<GenType,std::map<int,TGraph*>> map_map_gr_wire;
  std::map<GenType,TGraph*> map_h2_BackTracked;
  std::map<GenType,TH2D*> map_h2_PosHit;
  Helper h;
  TLatex* NHit           ;
  TLatex* NHitTot        ;
  TLatex* NeutrinoEnergy ;
  TLatex* VertexXPosition;

};


#endif
