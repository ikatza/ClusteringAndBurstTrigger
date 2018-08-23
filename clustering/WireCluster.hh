#ifndef WIRECLUSTER_HH
#define WIRECLUSTER_HH
#include "WireHit.hh"

class WireCluster {

public:
  WireCluster(std::vector<WireHit*>& cHitVector);
  WireCluster():
    fStartChannel  (0),
    fEndChannel    (0),
    fNChannel      (0),
    fChannelWidth  (0),
    fNHits         (0),
    fType          (0),
    fTriggerFlag   (0),
    fMarleyIndex   (0),
    fCutN          (0),
    fIsSelected    (0),
    fHitSADC       (0),
    fFirstHitTime  (0),
    fLastHitTime   (0),
    fTimeWidth     (0),
    fMC_UnderlyingE(0),
    fNElectron     (0),
    fEReco         (0),
    fHitVector     (),
    fPurity        (){
    fTruePosition[0] = 0;
    fTruePosition[1] = 0;
    fTruePosition[2] = 0;
    fRecoPosition[0] = 0;
    fRecoPosition[1] = 0;
    fRecoPosition[2] = 0;
  };
    
  ~WireCluster(){
    fHitVector.clear();
  };

  int    GetStartChannel  () const { return fStartChannel  ; };
  int    GetEndChannel    () const { return fEndChannel    ; };
  int    GetNChannel      () const { return fNChannel      ; };
  int    GetChannelWidth  () const { return fChannelWidth  ; };
  int    GetNHit          () const { return fNHits         ; };
  int    GetType          () const { return fType          ; };
  bool   GetTriggerFlag   () const { return fTriggerFlag   ; };
  int    GetFailedCut     () const { return fCutN          ; };
  bool   GetIsSelected    () const { return fIsSelected    ; };
  float  GetHitSumADC     () const { return fHitSADC       ; };
  float  GetFirstTimeHit  () const { return fFirstHitTime  ; };
  float  GetLastTimeHit   () const { return fLastHitTime   ; };
  float  GetTimeWidth     () const { return fTimeWidth     ; };
  int    GetMarleyIndex   () const { return fMarleyIndex   ; };
  double GetMC_UnderlyingE() const { return fMC_UnderlyingE; };
  double GetNElectron     () const { return fNElectron     ; };
  double GetEReco         () const { return fEReco         ; };
  double GetTruePosition(const int index=0) const { return fTruePosition[index]; };
  double GetRecoPosition(const int index=0) const { return fRecoPosition[index]; };
  double GetPurity      (GenType gen=kSNnu) const { return fPurity.at(gen); };
  std::vector<WireHit*> GetHitVector() const { return fHitVector; };

  void SetHitSADC     (const float  cHitSADC     ) { fHitSADC      = cHitSADC     ; };
  void SetTriggerFlag (const int    cTriggerFlag ) { fTriggerFlag  = cTriggerFlag ; };
  void SetIsSelected  (const bool   cIsSelected  ) { fIsSelected   = cIsSelected  ; };
  void SetFailedCut   (const int    cCutN        ) { fCutN         = cCutN        ; };
  void SetEReco       (const double cE           ) { fEReco        = cE           ; };
  void Print(const bool printHit=false) const;

  friend bool goodClusterHits(const WireCluster*, int);

private:
  int                      fStartChannel   ;
  int                      fEndChannel     ;
  int                      fNChannel       ;
  int                      fChannelWidth   ;
  int                      fNHits          ;
  int                      fType           ;
  bool                     fTriggerFlag    ;
  int                      fMarleyIndex    ;
  int                      fCutN           ;
  bool                     fIsSelected     ;
  float                    fHitSADC        ;
  float                    fFirstHitTime   ;
  float                    fLastHitTime    ;
  float                    fTimeWidth      ;
  double                   fMC_UnderlyingE ;
  double                   fNElectron      ;
  double                   fEReco          ;
  double                   fTruePosition[3];
  double                   fRecoPosition[3];
  std::vector<WireHit*>    fHitVector      ;
  std::map<GenType,double> fPurity         ;

public:
//  ClassDef(WireCluster,1)  //Simple event class

};


inline bool goodClusterHits(const WireCluster *c, int n){
  return c->fNHits>=n;
};

#endif
