#ifndef OPTICALHIT_HH
#define OPTICALHIT_HH

#include <iostream>

class OpticalHit {
public:
  OpticalHit(int GenType,
             double X, double Y, double Z,
             double PeakTime, double Width, double PE, int chan);
  OpticalHit():
    fGenType(0),
    fTime   (0),
    fWidth  (0),
    fPE     (0),
    fMarleyIndex(-1),
    fChannel(0)
    {
      fRecoPosition[0]=0;
      fRecoPosition[1]=0;
      fRecoPosition[2]=0;
    };

  ~OpticalHit(){};
  void   Print          () const {std::cout << "fGenType " << fGenType << std::endl; };
  double GetGenType     () const { return fGenType     ; };
  double GetTime        () const { return fTime        ; };
  double GetWidth       () const { return fWidth       ; };
  double GetChannel     () const { return fChannel     ; };
  double GetSPE         () const { return fPE          ; };
  double GetRecoPosition(const int i=0) const { return fRecoPosition[i]; };
  void   SetMarleyIndex (const int i=0) { fMarleyIndex = i;};
  int    GetMarleyIndex () const { return fMarleyIndex;};
  
  friend bool OpticalHitOrderedInTimePtr (const OpticalHit*, const OpticalHit*);
  friend bool OpticalHitOrderedInSpacePtr(const OpticalHit*, const OpticalHit*);
  friend bool OpticalHitOrderedInTime    (const OpticalHit&, const OpticalHit&);
  friend bool OpticalHitOrderedInSpace   (const OpticalHit&, const OpticalHit&);
  
private:
  int    fGenType;
  double fRecoPosition[3];
  double fTime;
  double fWidth;
  double fPE;
  int    fMarleyIndex;
  int    fChannel;
public:
  // ClassDef(OpticalHit,1)
};

inline bool OpticalHitOrderedInTimePtr (const OpticalHit* lhs, const OpticalHit* rhs){ return (lhs->fTime           <rhs->fTime);            };
inline bool OpticalHitOrderedInSpacePtr(const OpticalHit* lhs, const OpticalHit* rhs){ return (lhs->fRecoPosition[2]<rhs->fRecoPosition[2]); };
inline bool OpticalHitOrderedInTime    (const OpticalHit& lhs, const OpticalHit& rhs){ return (lhs.fTime            <rhs.fTime);             };
inline bool OpticalHitOrderedInSpace   (const OpticalHit& lhs, const OpticalHit& rhs){ return (lhs.fRecoPosition[2] <rhs.fRecoPosition[2]);  };


#endif
