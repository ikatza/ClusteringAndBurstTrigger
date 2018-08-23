#include "SmartTrigger.hh"

#include "Helper.h"

SmartTrigger::SmartTrigger(const std::string f):
  fLikelihood_Sign(NULL),
  fLikelihood_Back(NULL),
  fLookupDiscriminator(NULL),
  fFeature(f),
  fConfig(-1) {
  
};


void SmartTrigger::SetResultFromCache(const std::string f) {
  TFile* file = new TFile(f.c_str(), "READ");
  fLikelihood_Sign = (TH1D*)file->Get("Likelihood_Sign");
  fLikelihood_Back = (TH1D*)file->Get("Likelihood_Back");
  fLookupDiscriminator = (TH1D*)file->Get("LookupDiscriminator");

  fLikelihood_Sign    ->SetDirectory(NULL);
  fLikelihood_Back    ->SetDirectory(NULL);
  fLookupDiscriminator->SetDirectory(NULL);

  file->Close();
};


void SmartTrigger::SetBinning(const std::vector<double> binning) {

  int nbin = (int)binning.size()-1;
  if (nbin <= 0) {
    exit(0);
  }

  double* bin_d = new double[nbin+1];
  int i=0;
  for (auto const& b : binning) {
    bin_d[i++] = b;
  }

  fLikelihood_Sign = new TH1D("Likelihood_Sign", (";"+fFeature).c_str(),nbin, bin_d);
  fLikelihood_Back = new TH1D("Likelihood_Back", (";"+fFeature).c_str(),nbin, bin_d);
  fLookupDiscriminator = new TH1D("LookupDiscriminator", (";"+fFeature).c_str(),nbin, bin_d);
  fLikelihood_Sign->Sumw2();
  fLikelihood_Back->Sumw2();
  fLookupDiscriminator->Sumw2();
};


void SmartTrigger::ConstructLikelihood(const std::string InputFile,
                                       const std::string OutputFile) {
  
  if (fConfig == -1) {
    std::cout << "Configuration wasn't set" << std::endl;
    exit(0);
  }

  TFile InFile(InputFile.c_str(), "READ");
  TTree* Tree = (TTree*)InFile.Get("ClusteredWireHit");

  //TBranch** br = nullptr;
  double feat, type;
  int config;
  Tree->SetBranchAddress(fFeature.c_str(), &feat);
  Tree->SetBranchAddress("Config", &config);
  Tree->SetBranchAddress("Type",   &type);

  for (int i=0; i<Tree->GetEntries(); ++i) {
    PrintProgress(i,Tree->GetEntries());
    Tree->GetEntry(i);

    if (config != fConfig) continue;
    if (type == 1) { fLikelihood_Sign->Fill(feat); }
    else           { fLikelihood_Back->Fill(feat); }
  }

  fLikelihood_Sign->Scale(1./fLikelihood_Sign->Integral());
  fLikelihood_Back->Scale(1./fLikelihood_Back->Integral());

  for(int i=0; i<=fLikelihood_Sign->GetXaxis()->GetNbins(); ++i){
    double Back     = fLikelihood_Back->GetBinContent(i);
    double SignBack = fLikelihood_Sign->GetBinContent(i)+Back;
    if (Back > 0) {
      fLookupDiscriminator->SetBinContent(i, 2.*(Back + SignBack * (-1 + TMath::Log(SignBack / Back))));
    }
    else if (Back <= 0 && SignBack <= 0) {
      fLookupDiscriminator->SetBinContent(i, 0);
    }
    else if (Back <= 0) {
      fLookupDiscriminator->SetBinContent(i, 10000);
    }
  }
  std::cout << "Creating " << OutputFile << std::endl;
  TFile OutFile(OutputFile.c_str(), "RECREATE");
  fLikelihood_Sign->Write();
  fLikelihood_Back->Write();
  fLookupDiscriminator->Write();
  OutFile.Close();
  InFile .Close();

};



