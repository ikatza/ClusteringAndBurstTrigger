#include "OpticalHit.hh"

//ClassImp(OpticalHit)

OpticalHit::OpticalHit(int GenType,
                       double X, double Y, double Z,
                       double PeakTime, double Width, double PE, int chan):
  fGenType(GenType),
  fTime(PeakTime),
  fWidth(Width),
  fPE(PE),
  fChannel(chan)
{
  fRecoPosition[0] = (X);
  fRecoPosition[1] = (Y);
  fRecoPosition[2] = (Z);
  };
