#ifndef WIREPOSITIONTIMINGDISPLAY_HH
#define WIREPOSITIONTIMINGDISPLAY_HH

#include "Display.hh"

#include "TH2D.h"
#include "TPad.h"

template <class T>
class Boundary {
public:
  T low;
  T high;
  Boundary(const T l=0, const T h=1):low(l),high(h) {};
  T Spread() const { return high - low; };
  bool IsIn(const T t) const { return (low<t) && (t<high); }
  bool IsInOrOnBorder(const T t) const { return (low<=t) && (t<=high); }
};


class WirePositionTimingDisplay: public Display {

public:
  WirePositionTimingDisplay(const std::string, const std::string);
  void DisplayEvent(const int, const int);
  void LookForAPA(const GenType);

  int fAPA;
  int fChan;
  double fTime;
  std::map<GenType,TH2D*> f_map_wire_time;
  std::map<size_t,Boundary<int>> fAPA_Bounds;
};
  
#endif
