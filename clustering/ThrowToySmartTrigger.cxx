#include "TCanvas.h"
#include "TH1D.h"
#include "TPad.h"
#include "TGraph.h"


#include "SmartTriggerToyThrower.hh"


int main(int argc, char** argv) {

  int opt;
  std::string InputFile = "";
  std::string OutputFile = "";
  int nNeutrino=10, nToy=1000;
  double backRate=1;
  // Shut GetOpt error messages down (return '?'): 
  extern char *optarg;
  extern int  optopt;
  // Retrieve the options:
  while ( (opt = getopt(argc, argv, "t:c:f:i:o:b:n:")) != -1 ) {  // for each option...
    switch ( opt ) {
    case 'i':
      InputFile = optarg;
      break;
    case 'o':
      OutputFile = optarg;
      break;
    case 'n':
      nNeutrino = atof(optarg);
      break;
    case 'b':
      backRate = std::stod(optarg);
      break;
    case 't':
      nToy = atof(optarg);
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

  SmartTriggerToyThrower sttt(InputFile);
  sttt.SetNToys(nToy);
  sttt.SetBackgroundRateAndNumberOfNeutrino(backRate,nNeutrino);
  sttt.ThrowToys();

  TH1D* disc_sign = sttt.GetDiscriminator_Sign();
  TH1D* disc_back = sttt.GetDiscriminator_Back();
  TH1D* norm_sign = sttt.GetNumber_Sign();
  TH1D* norm_back = sttt.GetNumber_Back();

  norm_sign->SetLineWidth(2);
  norm_back->SetLineWidth(2);
  disc_sign->SetLineWidth(2);
  disc_back->SetLineWidth(2);

  norm_sign->SetLineColor(kRed);
  norm_back->SetLineColor(kBlue);
  disc_sign->SetLineColor(kRed);
  disc_back->SetLineColor(kBlue);

  double max = std::max(norm_sign->GetMaximum(),
                        norm_back->GetMaximum());
  norm_sign->SetMaximum(max*1.2);
  norm_back->SetMaximum(max*1.2);

  max = std::max(disc_sign->GetMaximum(),
                 disc_back->GetMaximum());
  disc_sign->SetMaximum(max*1.2);
  disc_back->SetMaximum(max*1.2);

  TCanvas c;
  c.Print("SmartTrigger.pdf[");

  TH1D* likelihood_Sign = sttt.GetSmartTrigger().GetLikelihood_Signal();
  TH1D* likelihood_Back = sttt.GetSmartTrigger().GetLikelihood_Background();
  likelihood_Sign->SetLineWidth(2);
  likelihood_Back->SetLineWidth(2);
  likelihood_Sign->SetLineColor(kRed);
  likelihood_Back->SetLineColor(kBlue);
  max = std::max(likelihood_Sign->GetMaximum(),
                 likelihood_Back->GetMaximum());
  likelihood_Sign->SetMaximum(max*1.2);
  likelihood_Back->SetMaximum(max*1.2);
  likelihood_Sign->SetStats(0);
  likelihood_Back->SetStats(0);
  likelihood_Sign->Draw();
  likelihood_Back->Draw("SAME");
  gPad->RedrawAxis();
  c.Print("SmartTrigger.pdf");

  TH1D* disc = sttt.GetSmartTrigger().GetLookupDiscriminator();
  disc->SetStats(0);
  disc->SetLineWidth(2);
  disc->SetLineColor(kRed);
  disc->Draw();
  gPad->SetLogy();
  gPad->RedrawAxis();
  c.Print("SmartTrigger.pdf");
    
  norm_sign->Draw();
  norm_back->Draw("SAME");
  gPad->RedrawAxis();
  c.Print("SmartTrigger.pdf");

  disc_sign->Draw();
  disc_back->Draw("SAME");
  gPad->RedrawAxis();
  c.Print("SmartTrigger.pdf");

  TGraph* ROC_Curve_norm = new TGraph(norm_back->GetXaxis()->GetNbins());
  TGraph* ROC_Curve_disc = new TGraph(disc_back->GetXaxis()->GetNbins());
  
  double norm_back_int = norm_back->Integral();
  double disc_back_int = disc_back->Integral();
  double norm_sign_int = norm_sign->Integral();
  double disc_sign_int = disc_sign->Integral();
  
  for (int i=0; i<=norm_back->GetXaxis()->GetNbins(); ++i) {
    double eff  = norm_sign->Integral(i,norm_back->GetXaxis()->GetNbins()+1)/norm_sign_int;
    double rate = backRate * norm_back->Integral(i,norm_back->GetXaxis()->GetNbins()+1)/norm_back_int;
    ROC_Curve_norm->SetPoint(i,eff,rate);
  }
  
  for (int i=0; i<=disc_back->GetXaxis()->GetNbins(); ++i) {
    double eff  = disc_sign->Integral(i,disc_back->GetXaxis()->GetNbins()+1)/disc_sign_int;
    double rate = backRate * disc_back->Integral(i,disc_back->GetXaxis()->GetNbins()+1)/disc_back_int;
    ROC_Curve_disc->SetPoint(i,eff,rate);
  }

  ROC_Curve_norm->SetLineColor  (kBlue);
  ROC_Curve_norm->SetMarkerColor(kBlue);
  ROC_Curve_disc->SetLineColor  (kRed);
  ROC_Curve_disc->SetMarkerColor(kRed);

  gPad->SetLogy(false);
  ROC_Curve_norm->SetTitle(";Efficiency;Background rate [Hz]");
  ROC_Curve_norm->Draw("ALP");
  ROC_Curve_disc->Draw("LP");
  gPad->RedrawAxis();
  c.Print("SmartTrigger.pdf");

  c.Print("SmartTrigger.pdf]");
}
