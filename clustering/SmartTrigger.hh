#ifndef SMARTTRIGGER_HH
#define SMARTTRIGGER_HH
#include <algorithm>
#include <iostream>

#include "TAxis.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1D.h"
#include "TTree.h"

#include "WireCluster.hh"

class SmartTrigger {
public:
  SmartTrigger(const std::string);
  ~SmartTrigger(){
    if (fLikelihood_Sign    ) delete fLikelihood_Sign    ;
    if (fLikelihood_Back    ) delete fLikelihood_Back    ;
    if (fLookupDiscriminator) delete fLookupDiscriminator;
    fLikelihood_Sign = NULL;
    fLikelihood_Back = NULL;
    fLookupDiscriminator = NULL;
  };

  void ConstructLikelihood(const std::string, const std::string);
  void SetBinning(const std::vector<double>);
  void SetResultFromCache(const std::string);

  void SetConfig(const int c=-1) { fConfig = c; };
  int  GetConfig() const { return fConfig; };
  
  TH1D* GetLikelihood_Signal    () const { return fLikelihood_Sign; };
  TH1D* GetLikelihood_Background() const { return fLikelihood_Back; };
  TH1D* GetLookupDiscriminator  () const { return fLookupDiscriminator; };
  
  double GetEstimator(const WireCluster& cluster) const {
    return cluster.GetEReco();
  };
  

  double GetDiscriminator(const double feat) const {
    if (!fLookupDiscriminator) {
      std::cout << "Likelihoods were not instantiated!" << std::endl;
      exit(1);
    }
    int bin = fLookupDiscriminator->FindBin(feat);
    bin = std::min(bin,fLookupDiscriminator->GetXaxis()->GetNbins());
    return fLookupDiscriminator->GetBinContent(bin);
  };
  
  // double GetDiscriminatorVecDouble(const std::vector<double> feat) const {
  //   double r=0;
  //   for (std::vector<double>iterator it = feat.begin(); it != feat.end(); ++it) {
  //     r += GetDiscriminator(it);
  //   }
  //   return r;
  // };

  double GetDiscriminator(const WireCluster& cluster) const {
    return GetDiscriminator(GetEstimator(cluster));
  };

  // double GetDiscriminatorVecCluster(const std::vector<WireCluster&> cluster) const {
  //   double r = 0;
  //   for (auto const& it : cluster) {
  //     r += GetDiscriminator(GetEstimator(it));
  //   }
  //   return r;
  // };
  
private:
  TH1D* fLikelihood_Sign;
  TH1D* fLikelihood_Back;
  TH1D* fLookupDiscriminator;
  std::string fFeature;
  int fConfig; 
  
};

#endif
