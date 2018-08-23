#include "SmartTriggerToyThrower.hh"


void SmartTriggerToyThrower::ThrowToys() {

  std::vector<double> discriminator_Sign_vec;
  std::vector<double> discriminator_Back_vec;
  std::vector<double> number_Sign_vec;
  std::vector<double> number_Back_vec;
  if (fMeanBack<=0 || fMeanSign<=0) {
    std::cout << "You should call SetBackgroundRateAndNumberOfNeutrino before starting to throw toys"<< std::endl;
    exit(1);
  }
  for (size_t toy=0; toy<fNToys; ++toy) {
    PrintProgress(toy,fNToys);
    double norm_back  = gRandom->PoissonD(fMeanBack);
    double norm_back2 = gRandom->PoissonD(fMeanBack);
    double norm_sign  = gRandom->PoissonD(fMeanSign);
    number_Back_vec.push_back(norm_back);
    number_Sign_vec.push_back(norm_back2+norm_sign);
    double feat_back = 0;
    double feat_sign = 0;
    for (int i=0; i<norm_back; ++i) {
      feat_back += fSmTrig.GetDiscriminator(fLikelihood_Back->GetRandom());
    }
    for (int i=0; i<norm_sign; ++i) {
      feat_sign += fSmTrig.GetDiscriminator(fLikelihood_Sign->GetRandom());
    }
    for (int i=0; i<norm_back2; ++i) {
      feat_sign += fSmTrig.GetDiscriminator(fLikelihood_Back->GetRandom());
    }
    discriminator_Sign_vec.push_back(feat_sign);
    discriminator_Back_vec.push_back(feat_back);
  }

  double max_b = *max_element(discriminator_Back_vec.begin(), discriminator_Back_vec.end());
  double min_b = *min_element(discriminator_Back_vec.begin(), discriminator_Back_vec.end());
  double max_s = *max_element(discriminator_Sign_vec.begin(), discriminator_Sign_vec.end());
  double min_s = *min_element(discriminator_Sign_vec.begin(), discriminator_Sign_vec.end());
  double max = max_b;
  double min = min_b;//std::min(0, min_b);
  fDiscriminator_Sign = new TH1D("Discriminator_Sign", ";Discriminator;nEvents", 1000, min, max*2.);
  fDiscriminator_Back = new TH1D("Discriminator_Back", ";Discriminator;nEvents", 1000, min, max*2.);
  for (auto const& it: discriminator_Sign_vec) fDiscriminator_Sign->Fill(it);
  for (auto const& it: discriminator_Back_vec) fDiscriminator_Back->Fill(it);
  AddOverflow(fDiscriminator_Sign);
  AddOverflow(fDiscriminator_Back);
  fDiscriminator_Sign->SetStats(0);
  fDiscriminator_Back->SetStats(0);
  
  max_b = *max_element(number_Back_vec.begin(), number_Back_vec.end());
  min_b = *min_element(number_Back_vec.begin(), number_Back_vec.end());
  max_s = *max_element(number_Sign_vec.begin(), number_Sign_vec.end());
  min_s = *min_element(number_Sign_vec.begin(), number_Sign_vec.end());
  max = max_b*2.;
  min = min_b;//std::min(0, min_b);
  fNumber_Sign = new TH1D("Number_Sign", ";Number of clusters;nEvents", max-min, min, max);
  fNumber_Back = new TH1D("Number_Back", ";Number of clusters;nEvents", max-min, min, max);
  for (auto const& it: number_Sign_vec) fNumber_Sign->Fill(it);
  for (auto const& it: number_Back_vec) fNumber_Back->Fill(it);
  AddOverflow(fNumber_Sign);
  AddOverflow(fNumber_Back);
  fNumber_Sign->SetStats(0);
  fNumber_Back->SetStats(0);

  fToyThrown = true;

};
