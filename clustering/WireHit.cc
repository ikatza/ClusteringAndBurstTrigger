#include "WireHit.hh"

//ClassImp(WireHit)

void WireHit::Print() const {
  std::cout << "fHitView " << fHitView
            << ",  fGenType " << fGenType
            << ",  fHitChan " << fHitChan
            << ",  fHitSADC " << fHitSADC
            << ",  fHitRMS  " << fHitRMS 
            << ",  fHitTime " << fHitTime
            << ",  fTrDepositedE " << fTrDepositedE
            << ",  fTrEPart      " << fTrEPart     
            << ",  fTrIDPart     " << fTrIDPart << std::endl;
  // std::cout << "fTruePosition[0] " << fTruePosition[0] << std::endl;
  // std::cout << "fTruePosition[1] " << fTruePosition[1] << std::endl;
  // std::cout << "fTruePosition[2] " << fTruePosition[2] << std::endl;
  // std::cout << "fRecoPosition[0] " << fRecoPosition[0] << std::endl;
  // std::cout << "fRecoPosition[1] " << fRecoPosition[1] << std::endl;
  // std::cout << "fRecoPosition[2] " << fRecoPosition[2] << std::endl;
};
