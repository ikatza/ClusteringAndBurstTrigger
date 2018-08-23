#ifndef OPTICALCLUSTER_HH
#define OPTICALCLUSTER_HH
#include <iostream>
#include <map>
#include <vector>

#include "OpticalHit.hh"
#include "WireHit.hh"

class OpticalCluster {
public:
  OpticalCluster(std::vector<OpticalHit*>& cHitVector);
  OpticalCluster():
    fEvent         (0),
    fNChannel      (0),
    fNHit          (0),
    fType          (0),
    fHitSPE        (0),
    fMarleyIndex   (0),
    fFirstHitTime  (0),
    fLastHitTime   (0),
    fTimeWidth     (0),
    fEnu           (0),
    fHitVector     (),
    fPurity        (){
    fRecoPosition[0] = 0;
    fRecoPosition[1] = 0;
    fRecoPosition[2] = 0;
    fRecoWidth   [0] = 0;
    fRecoWidth   [1] = 0;
    fRecoWidth   [2] = 0;
    fTrueVertexPos[0] = 0;
    fTrueVertexPos[1] = 0;
    fTrueVertexPos[2] = 0;
    fTrueVertexPos[3] = 0;
    fTrueVertexDir[0] = 0;
    fTrueVertexDir[1] = 0;
    fTrueVertexDir[2] = 0;
  };
  ~OpticalCluster(){
    fHitVector.clear();
    fPurity.clear();
  };

  std::vector<WireHit*>* CreateFakeWireHits() {
    return new std::vector<WireHit*>();
  };

  int                      GetEvent       ()              const { return fEvent           ; };
  int                      GetNChannel    ()              const { return fNChannel        ; };
  int                      GetNHit        ()              const { return fNHit            ; };
  int                      GetType        ()              const { return fType            ; };
  int                      GetHitSPE      ()              const { return fHitSPE          ; };
  int                      GetMarleyIndex ()              const { return fMarleyIndex     ; };
  float                    GetFirstHitTime()              const { return fFirstHitTime    ; };
  float                    GetLastHitTime ()              const { return fLastHitTime     ; };
  float                    GetTimeWidth   ()              const { return fTimeWidth       ; };
  float                    GetRecoPosition(const int i=0) const { return fRecoPosition[i] ; };
  float                    GetRecoWidth   (const int i=0) const { return fRecoWidth   [i] ; };
  double                   GetTrueEnu     ()              const { return fEnu             ; };
  double                   GetTrueVertPos (const int i=0) const { return fTrueVertexPos[i]; };
  double                   GetTrueVertDir (const int i=0) const { return fTrueVertexDir[i]; };
  std::vector<OpticalHit*> GetHitVector   ()              const { return fHitVector       ; };

  void Print() const;
  friend bool OpticalClusterOrderedInTimePtr (const OpticalCluster*, const OpticalCluster*);
private:
  int    fEvent           ;
  int    fNChannel        ;
  int    fNHit            ;
  int    fType            ;
  int    fHitSPE          ;
  int    fMarleyIndex     ;
  float  fFirstHitTime    ;
  float  fLastHitTime     ;
  float  fTimeWidth       ;
  float  fRecoPosition [3];
  float  fRecoWidth    [3];
  double fEnu             ;
  double fTrueVertexPos[4];
  double fTrueVertexDir[3];
  std::vector<OpticalHit*> fHitVector;
  std::map<GenType,double> fPurity;

public:
  // ClassDef(OpticalCluster,1)

};

inline bool OpticalClusterOrderedInTimePtr (const OpticalCluster* lhs, const OpticalCluster* rhs){ return (lhs->fFirstHitTime<rhs->fFirstHitTime); };

#endif
