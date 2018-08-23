#ifndef WIREHIT_HH
#define WIREHIT_HH
#include "Helper.h"

class WireHit {
public:
  WireHit(int   cHitView,      int   cGenType, int   cHitChan,
          float cHitTime,      float cHitSADC, float cHitRMS,
          float cTrDepositedE, float cTrEPart, int   cTrIDPart,
          float cRecPosX,      float cRecPosY, float cRecPosZ,
          float cTrPosX,       float cTrPosY,  float cTrPosZ,
          int   cMarleyIndex,  float cNelec):
    fHitView(cHitView),
    fGenType(cGenType),
    fHitChan(cHitChan),
    fHitSADC(cHitSADC),
    fHitRMS (cHitRMS ),
    fHitTime(cHitTime),
    fMarleyIndex(cMarleyIndex),
    fTrDepositedE(cTrDepositedE),
    fNElectron(cNelec),
    fTrEPart(cTrEPart),
    fTrIDPart(cTrIDPart),
    fHitChanNorm(0),
    fHitTimeNorm(0){
    fTruePosition[0] = cTrPosX;
    fTruePosition[1] = cTrPosY;
    fTruePosition[2] = cTrPosZ;
    fRecoPosition[0] = cRecPosX;
    fRecoPosition[1] = cRecPosY;
    fRecoPosition[2] = cRecPosZ;
  };
  
  WireHit():
    fHitView(0),
    fGenType(0),
    fHitChan(0),
    fHitSADC(0),
    fHitRMS (0),
    fHitTime(0),
    fMarleyIndex(0),
    fTrDepositedE(0),
    fNElectron(0),
    fTrEPart (0),
    fTrIDPart(0),
    fHitChanNorm(0),
    fHitTimeNorm(0){
    fTruePosition[0] = (0);
    fTruePosition[1] = (0);
    fTruePosition[2] = (0);
    fRecoPosition[0] = (0);
    fRecoPosition[1] = (0);
    fRecoPosition[2] = (0);
  };
  
  WireHit(const WireHit&hit):
    fHitView(hit.fHitView),
    fGenType(hit.fGenType),
    fHitChan(hit.fHitChan),
    fHitSADC(hit.fHitSADC),
    fHitRMS (hit.fHitRMS ),
    fHitTime(hit.fHitTime),
    fMarleyIndex (hit.fMarleyIndex),
    fTrDepositedE(hit.fTrDepositedE),
    fNElectron(hit.fNElectron),
    fTrEPart     (hit.fTrEPart ),
    fTrIDPart    (hit.fTrIDPart),
    fHitChanNorm(hit.fHitChanNorm),
    fHitTimeNorm(hit.fHitTimeNorm){
    fTruePosition[0] = hit.fTruePosition[0];
    fTruePosition[1] = hit.fTruePosition[1];
    fTruePosition[2] = hit.fTruePosition[2];
    fRecoPosition[0] = hit.fRecoPosition[0];
    fRecoPosition[1] = hit.fRecoPosition[1];
    fRecoPosition[2] = hit.fRecoPosition[2];
  };

  ~WireHit(){};
    
  friend bool WireHitOrderedInTime    (const WireHit& lhs, const WireHit& rhs);
  friend bool WireHitOrderedInSpace   (const WireHit& lhs, const WireHit& rhs);
  friend bool WireHitOrderedInTimePtr (const WireHit* lhs, const WireHit* rhs);
  friend bool WireHitOrderedInSpacePtr(const WireHit* lhs, const WireHit* rhs);

  friend double GetMeanChannel  (const std::vector<WireHit*> hits);
  friend double GetStdDevChannel(const std::vector<WireHit*> hits);
  friend double GetMeanTime     (const std::vector<WireHit*> hits);
  friend double GetStdDevTime   (const std::vector<WireHit*> hits);
  
  friend void NormaliseChannelAndTime(std::vector<WireHit*> hits);

//  void  setOptical(const bool op=true) { fIsOptical=op; };
  
  int   GetHitView     () const { return fHitView;      };
  int   GetGenType     () const { return fGenType;      };
  int   GetHitChan     () const { return fHitChan;      };
  int   GetMarleyIndex () const { return fMarleyIndex;  };
  float GetHitTime     () const { return fHitTime;      };
  float GetHitSADC     () const { return fHitSADC;      };
  float GetHitRMS      () const { return fHitRMS;       };
  float GetTrDepositedE() const { return fTrDepositedE; };
  float GetNElectron   () const { return fNElectron;    };
  float GetTrEPart     () const { return fTrEPart;      };
  float GetTrIDPart    () const { return fTrIDPart;     };
  float GetHitChanNorm () const { return fHitChanNorm;  };
  float GetHitTimeNorm () const { return fHitTimeNorm;  };
  float GetRecoPosition(const int index=0) const { return fRecoPosition[index]; };
  float GetTruePosition(const int index=0) const { return fTruePosition[index]; };
  
  void Print() const;

  
private:
  int   fHitView;
  int   fGenType;
  int   fHitChan;
  float fHitSADC;
  float fHitRMS ;
  float fHitTime;
  float fRecoPosition[3];
  int   fMarleyIndex ;
  float fTrDepositedE;
  float fNElectron   ;
  float fTrEPart     ;
  int   fTrIDPart    ;
  float fTruePosition[3];
  
  mutable float fHitChanNorm;
  mutable float fHitTimeNorm;

public:
  //ClassDef(WireHit,1)  //Simple event class
};

inline bool WireHitOrderedInTimePtr (const WireHit* lhs, const WireHit* rhs){ return (lhs->fHitTime<rhs->fHitTime); };
inline bool WireHitOrderedInSpacePtr(const WireHit* lhs, const WireHit* rhs){ return (lhs->fHitChan<rhs->fHitChan); };
inline bool WireHitOrderedInTime    (const WireHit& lhs, const WireHit& rhs){ return (lhs.fHitTime <rhs.fHitTime);  };
inline bool WireHitOrderedInSpace   (const WireHit& lhs, const WireHit& rhs){ return (lhs.fHitChan <rhs.fHitChan);  };

#endif
