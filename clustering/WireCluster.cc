
#include "WireCluster.hh"
#include "WireHit.hh"

//ClassImp(WireCluster)

//......................................................
WireCluster::WireCluster(std::vector<WireHit*>& vec_Hits):
  fStartChannel  (0),
  fEndChannel    (0),
  fNChannel      (0),
  fChannelWidth  (0),
  fNHits         (0),
  fType          (0),
  fTriggerFlag   (0),
  fMarleyIndex   (-1),
  fCutN          (0),
  fIsSelected    (0),
  fHitSADC       (0),
  fFirstHitTime  (0),
  fLastHitTime   (0),
  fTimeWidth     (0),
  fMC_UnderlyingE(0),
  fNElectron     (0),
  fHitVector     (),
  fPurity        (){
  
  fTruePosition[0] = 0;
  fTruePosition[1] = 0;
  fTruePosition[2] = 0;
  std::set<int> channels;
  fHitVector.clear();
  fStartChannel = 50000;
  fEndChannel   = 0;
  fFirstHitTime = 50000;
  fLastHitTime  = 0;
  std::map<int,int> map_index_number;

  fPurity[kOther   ] = 0;
  fPurity[kSNnu    ] = 0;
  fPurity[kAPA     ] = 0;
  fPurity[kCPA     ] = 0;
  fPurity[kAr39    ] = 0;
  fPurity[kNeutron ] = 0;
  fPurity[kKrypton ] = 0;
  fPurity[kPolonium] = 0;
  fPurity[kRadon   ] = 0;
  fPurity[kAr42    ] = 0;
  
  std::map<int,double>map_id_energy;
  for(size_t i = 0; i < vec_Hits.size(); ++i)
  {
    fHitVector.push_back(vec_Hits[i]);
    fHitSADC  += vec_Hits[i]->GetHitSADC(); 
    int   chan = vec_Hits[i]->GetHitChan();
    float time = vec_Hits[i]->GetHitTime(); 
    fNElectron += vec_Hits[i]->GetNElectron();
    channels.insert(chan);
    map_index_number[vec_Hits[i]->GetMarleyIndex()]++;
    if(chan<fStartChannel) fStartChannel=chan;
    if(chan>fEndChannel)   fEndChannel=chan;
    if(time<fFirstHitTime) fFirstHitTime=time;
    if(time>fLastHitTime)  fLastHitTime=time;

    if(fHitVector[i]->GetGenType()==1) fType++;
    fPurity[ConvertIntToGenType(fHitVector[i]->GetGenType())] += vec_Hits[i]->GetHitSADC();

    map_id_energy[vec_Hits[i]->GetTrIDPart()] = vec_Hits[i]->GetTrEPart()*1000.;
    fRecoPosition[0] += vec_Hits[i]->GetRecoPosition(0) * vec_Hits[i]->GetHitSADC();
    fRecoPosition[1] += vec_Hits[i]->GetRecoPosition(1) * vec_Hits[i]->GetHitSADC();
    fRecoPosition[2] += vec_Hits[i]->GetRecoPosition(2) * vec_Hits[i]->GetHitSADC();
    fTruePosition[0] += vec_Hits[i]->GetTruePosition(0) * vec_Hits[i]->GetHitSADC();
    fTruePosition[1] += vec_Hits[i]->GetTruePosition(1) * vec_Hits[i]->GetHitSADC();
    fTruePosition[2] += vec_Hits[i]->GetTruePosition(2) * vec_Hits[i]->GetHitSADC();
  }
  
  if(fHitSADC!=0) {
    fRecoPosition[0] /= fHitSADC;
    fRecoPosition[1] /= fHitSADC;
    fRecoPosition[2] /= fHitSADC;
    fTruePosition[0] /= fHitSADC;
    fTruePosition[1] /= fHitSADC;
    fTruePosition[2] /= fHitSADC;
  } else {
    fRecoPosition[0] /= -999;
    fRecoPosition[1] /= -999;
    fRecoPosition[2] /= -999;
    fTruePosition[0] /= -999;
    fTruePosition[1] /= -999;
    fTruePosition[2] /= -999;
  }

  
  fNHits = fHitVector.size();
  fNChannel = channels.size();

  //CALL THE CLUSTER MARLEY IF THERE ARE MORE THAN TWO MARLEY HITS IN IT.
  if (fType>=2) {
    fType = 1;

    int number_hits=0;
    for (auto const& it : map_index_number)
      if (it.second>number_hits) {
        number_hits = it.second;
        fMarleyIndex = it.first;
      }
  } else {
    fType = 0;
  }

  for (auto& it:fPurity) {
    // std::cout << GenName[it.first] << " " << it.second;
    it.second /= fHitSADC;
    // std::cout << "   -> " << it.second << std::endl;
  }
  fMC_UnderlyingE = 0;
  for (auto& it: map_id_energy) {
    fMC_UnderlyingE += it.second;
  }
  fChannelWidth = fEndChannel  - fStartChannel;
  fTimeWidth    = fLastHitTime - fFirstHitTime;

  //ORDER IN THE HIT VECTOR BY TIME.
  std::sort(fHitVector.begin(), fHitVector.end(),WireHitOrderedInTimePtr);
  
};

void WireCluster::Print(const bool printHit) const {
  std::cout << "*************************************" << std::endl;
  std::cout << "StartChan        " << fStartChannel    << std::endl;
  std::cout << "EndChan          " << fEndChannel      << std::endl;
  std::cout << "NChan            " << fNChannel        << std::endl;
  std::cout << "ChanWidth        " << fChannelWidth    << std::endl;
  std::cout << "NHits            " << fNHits           << std::endl;
  std::cout << "Type             " << fType            << std::endl;
  std::cout << "TriggerFlag      " << fTriggerFlag     << std::endl;
  std::cout << "HitSADC          " << fHitSADC         << std::endl;
  std::cout << "FirstHitTime     " << fFirstHitTime    << std::endl;
  std::cout << "LastHitTime      " << fLastHitTime     << std::endl;
  std::cout << "TimeWidth        " << fTimeWidth       << std::endl;
  std::cout << "MC_UnderlyingE   " << fMC_UnderlyingE  << std::endl;
  std::cout << "fTruePosition[0] " << fTruePosition[0] << std::endl;
  std::cout << "fTruePosition[1] " << fTruePosition[1] << std::endl;
  std::cout << "fTruePosition[2] " << fTruePosition[2] << std::endl;
  std::cout << "fRecoPosition[0] " << fRecoPosition[0] << std::endl;
  std::cout << "fRecoPosition[1] " << fRecoPosition[1] << std::endl;
  std::cout << "fRecoPosition[2] " << fRecoPosition[2] << std::endl;
  if(printHit)
    for(int i = 0; i < fNHits; i++)
      fHitVector[i]->Print();
  std::cout << "*************************************" << std::endl;
};
