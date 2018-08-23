#ifndef HITSIMPLEPLOT_HH
#define HITSIMPLEPLOT_HH
#include "InputManager.h"
#include "Helper.h"

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"


class HitSimplePlot {

public:
  HitSimplePlot();
  ~HitSimplePlot();
  void PlotRMS    (const int);
  void PlotWidth  (const int);
  void PlotRMS1D  (const int);
  void PlotWidth1D(const int);
  void PlotRMS2D  (const int, const int);
  void PlotWidth2D(const int, const int);

private:
  InputManager im;
  TCanvas* c;
  void   GetDimension(const int, int&, int&) const;
  double GetMean  (std::vector<double> const) const;
  double GetStdDev(std::vector<double> const) const;
  void   GetMinMax(std::vector<double> const, double&, double&) const;
  std::map<GenType, std::map<int,std::vector<double>>> GetMap(std::vector<double> const*,
                                                              std::vector<int> const*,
                                                              std::vector<int> const*) const;
  std::map<GenType, std::map<int,std::pair<std::vector<double>,std::vector<double>>>> Get2DMap(std::vector<double> const*,
                                                                                               std::vector<double> const*,
                                                                                               std::vector<int> const*,
                                                                                               std::vector<int> const*) const;
  void FillHistosRMS  (std::vector<double>*, std::map<GenType, TH1D*>&);
  void FillHistosWidth(std::vector<double>*, std::map<GenType, TH1D*>&);
  void Fill2DHistosRMS  (std::vector<double>*, std::vector<double>*, std::map<GenType, TH2D*>&);
  void Fill2DHistosWidth(std::vector<double>*, std::vector<double>*, std::map<GenType, TH2D*>&);

  std::map<int, int> nbin_1D_RMS;
  std::map<int, int> xmin_1D_RMS;
  std::map<int, int> xmax_1D_RMS;
  std::map<int, int> nbin_1D_Width;
  std::map<int, int> xmin_1D_Width;
  std::map<int, int> xmax_1D_Width;
  std::map<int, int> nbin_2D_x_RMS;
  std::map<int, int> xmin_2D_x_RMS;
  std::map<int, int> xmax_2D_x_RMS;
  std::map<int, int> nbin_2D_x_Width;
  std::map<int, int> xmin_2D_x_Width;
  std::map<int, int> xmax_2D_x_Width;
  std::map<int, int> nbin_2D_y_RMS;
  std::map<int, int> xmin_2D_y_RMS;
  std::map<int, int> xmax_2D_y_RMS;
  std::map<int, int> nbin_2D_y_Width;
  std::map<int, int> xmin_2D_y_Width;
  std::map<int, int> xmax_2D_y_Width;

};
#endif
