#ifndef SMARTTRIGGERTOYTHROWER_HH
#define SMARTTRIGGERTOYTHROWER_HH
#include "SmartTrigger.hh"

#include "Helper.h"

#include "TRandom3.h"

class SmartTriggerToyThrower {
public:
  SmartTriggerToyThrower(const std::string f):
    fSmTrig(""),
    fToyThrown(false),
    fNToys(0),
    fSeed(5056),
    fTimeWindow(10),
    fMeanBack(-1),
    fMeanSign(-1),
    fNumber_Sign(NULL),
    fNumber_Back(NULL),
    fDiscriminator_Sign(NULL),
    fDiscriminator_Back(NULL),
    fLikelihood_Sign(NULL),
    fLikelihood_Back(NULL) {
    
    delete gRandom;
    gRandom = new TRandom3(fSeed);

    fSmTrig.SetResultFromCache(f);
    fLikelihood_Sign = fSmTrig.GetLikelihood_Signal();
    fLikelihood_Back = fSmTrig.GetLikelihood_Background();
  }

  ~SmartTriggerToyThrower() {
    if (fNumber_Sign       ) delete fNumber_Sign;
    if (fNumber_Back       ) delete fNumber_Back;
    if (fDiscriminator_Sign) delete fDiscriminator_Sign;
    if (fDiscriminator_Back) delete fDiscriminator_Back;
    if (fLikelihood_Sign   ) delete fLikelihood_Sign;
    if (fLikelihood_Back   ) delete fLikelihood_Back;   
    fNumber_Sign = NULL;
    fNumber_Back = NULL;
    fLikelihood_Sign = NULL;
    fLikelihood_Back = NULL;
    fDiscriminator_Sign = NULL;
    fDiscriminator_Back = NULL;
  };
  
  void SetBackgroundRateAndNumberOfNeutrino(double backRate,
                                            int nNeutrino) {
    fMeanBack = backRate * fTimeWindow;
    fMeanSign = nNeutrino;
  };

  void ThrowToys();
  void   SetNToys(const size_t t=10000) { fNToys = t; };
  size_t GetNToys() const { return fNToys; };
  void   SetSeed(size_t s) {
    fSeed = s;
    delete gRandom;
    gRandom = new TRandom3(fSeed);
  };
  void SetTimeWindow(const double d=10) { fTimeWindow = d; };
  double GetTimeWindow() const { return fTimeWindow; };
  TH1D* GetDiscriminator_Sign() const { return fDiscriminator_Sign; };
  TH1D* GetDiscriminator_Back() const { return fDiscriminator_Back; };
  TH1D* GetNumber_Sign()        const { return fNumber_Sign; };
  TH1D* GetNumber_Back()        const { return fNumber_Back; };
  
  SmartTrigger& GetSmartTrigger() { return fSmTrig; };

private:
  SmartTrigger fSmTrig;
  bool fToyThrown;
  size_t fNToys;
  size_t fSeed;
  double fTimeWindow;
  double fMeanBack;
  double fMeanSign;
  TH1D* fNumber_Sign;
  TH1D* fNumber_Back;
  TH1D* fDiscriminator_Sign;
  TH1D* fDiscriminator_Back;
  TH1D* fLikelihood_Sign;
  TH1D* fLikelihood_Back;

};

#endif
