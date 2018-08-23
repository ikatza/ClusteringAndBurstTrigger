#include <unistd.h>
#include <iostream>

#include "TH1D.h"
#include "TCanvas.h"

#include "SmartTrigger.hh"



int main(int argc, char** argv) {

  int opt;
  int Config=-1;
  std::string InputFile = "";
  std::string OutputFile = "";
  std::string Feature = "";
// Shut GetOpt error messages down (return '?'): 
  extern char *optarg;
  extern int  optopt;
  // Retrieve the options:
  while ( (opt = getopt(argc, argv, "c:f:i:o:")) != -1 ) {  // for each option...
    switch ( opt ) {
    case 'c':
      Config = atoi(optarg);
      break;
    case 'i':
      InputFile = optarg;
      break;
    case 'o':
      OutputFile = optarg;
      break;
    case 'f':
      Feature = optarg;
      break;
    case '?':  // unknown option...
      std::cerr << "Unknown option: '" << char(optopt) << "'!" << std::endl;
      break;
    }
  }
  
  if (InputFile == "") {
    std::cout << "No input file!!" << std::endl;
    exit(1);
  }
 
  if (OutputFile == "") {
    std::cout << "No output file!!" << std::endl;
    exit(1);
  }
  
  std::cout << "Using file " << InputFile << std::endl;
  
  if (Config<0 || Config>5) {
    std::cout << "Invalid config " << Config << std::endl;
    exit(1);
  }

  SmartTrigger strig(Feature);
  strig.SetConfig(Config);
  std::map<std::string,std::vector<double>> bin;
  bin["SumADC"] = {0,1000,2000,3000,4000,5000,10000};
  bin["EReco"]  = {0,17,20,30,50};
  bin["NElectron"]  = {0,5000,10000,15000,20000,25000,30000,35000,40000,45000,70000,100000};
  if (bin.find(Feature) == bin.end()) {
    std::cout << "The requested feature dpesn't have a binning set, go to RunSmartTrigger.cxx and add one to the bin map." << std::endl;
    exit(1);
  }
  strig.SetBinning(bin.at(Feature));
  strig.ConstructLikelihood(InputFile, OutputFile);

  TH1D* Likelihood_Sign = strig.GetLikelihood_Signal    ();
  TH1D* Likelihood_Back = strig.GetLikelihood_Background();

  TCanvas c;
  AddOverflow(Likelihood_Sign);
  AddOverflow(Likelihood_Back);
  std::cout<< Likelihood_Back->GetBinContent(1)<< std::endl;
  Likelihood_Sign->SetMinimum(0);
  Likelihood_Back->SetMinimum(0);
  std::vector<double> max = {Likelihood_Sign->GetBinContent(Likelihood_Sign->GetMaximumBin()),
                             Likelihood_Back->GetBinContent(Likelihood_Back->GetMaximumBin())};
  double max_d = *max_element(max.begin(), max.end())*1.2;
  Likelihood_Sign->SetMaximum(max_d);
  Likelihood_Back->SetMaximum(max_d);
  Likelihood_Sign->SetLineColor(kRed);
  Likelihood_Back->SetLineColor(kBlue);
  Likelihood_Sign->Draw();
  Likelihood_Back->Draw("SAME");
  c.SaveAs("likelihoods.pdf");

  TH1D* Lookup = strig.GetLookupDiscriminator();
  Lookup->Draw();
  gPad->SetLogy();
  c.SaveAs("Lookup.pdf");



}



