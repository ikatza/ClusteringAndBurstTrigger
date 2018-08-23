#include "OpticalCluster.hh"



OpticalCluster::OpticalCluster(std::vector<OpticalHit*>& vec_Hits):
  fNChannel      (0),
  fNHit          (0),
  fType          (0),
  fHitSPE        (0),
  fFirstHitTime  (100000),
  fLastHitTime   (-100000),
  fTimeWidth     (0),
  fHitVector     (),
  fPurity        (){
  fRecoPosition[0] = 0;
  fRecoPosition[1] = 0;
  fRecoPosition[2] = 0;
  fRecoWidth   [0] = 0;
  fRecoWidth   [1] = 0;
  fRecoWidth   [2] = 0;

  double minPos[3] = { 10000,  10000,  10000};
  double maxPos[3] = {-10000, -10000, -10000};
  ///std::cout << vec_Hits.size() << std::endl;
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
  std::set<int> channels;
  fHitVector.clear();
  std::map<int,int> map_marleyevent_number;
  for(size_t i = 0; i < vec_Hits.size(); ++i) {
    fHitVector.push_back(vec_Hits[i]);
    fHitSPE   += vec_Hits[i]->GetSPE ();
    int   chan = vec_Hits[i]->GetChannel();
    float time = vec_Hits[i]->GetTime();
    map_marleyevent_number[vec_Hits[i]->GetMarleyIndex()]++;
    channels.insert(chan);

    // if(chan<fStartChannel) fStartChannel=chan;
    // if(chan>fEndChannel)   fEndChannel=chan;
    if(time<fFirstHitTime) fFirstHitTime=time;
    if(time>fLastHitTime ) fLastHitTime =time;

    if(fHitVector[i]->GetGenType()==1) fType++;

    fRecoPosition[0] += vec_Hits[i]->GetRecoPosition(0) * vec_Hits[i]->GetSPE();
    fRecoPosition[1] += vec_Hits[i]->GetRecoPosition(1) * vec_Hits[i]->GetSPE();
    fRecoPosition[2] += vec_Hits[i]->GetRecoPosition(2) * vec_Hits[i]->GetSPE();
    fPurity[ConvertIntToGenType(fHitVector[i]->GetGenType())] ++;//= vec_Hits[i]->getSPE();

    if(vec_Hits[i]->GetRecoPosition(0) < minPos[0]) minPos[0] = vec_Hits[i]->GetRecoPosition(0);
    if(vec_Hits[i]->GetRecoPosition(1) < minPos[1]) minPos[1] = vec_Hits[i]->GetRecoPosition(1);
    if(vec_Hits[i]->GetRecoPosition(2) < minPos[2]) minPos[2] = vec_Hits[i]->GetRecoPosition(2);
    if(vec_Hits[i]->GetRecoPosition(0) > maxPos[0]) maxPos[0] = vec_Hits[i]->GetRecoPosition(0);
    if(vec_Hits[i]->GetRecoPosition(1) > maxPos[1]) maxPos[1] = vec_Hits[i]->GetRecoPosition(1);
    if(vec_Hits[i]->GetRecoPosition(2) > maxPos[2]) maxPos[2] = vec_Hits[i]->GetRecoPosition(2);
  }
  int max=0;
  fMarleyIndex=0;
  for (auto const & it : map_marleyevent_number) {
    if (it.second>max) {
      fMarleyIndex=it.first;
      max = it.second;
    }
  }
  fRecoWidth[0] = std::abs(minPos[0] - maxPos[0]);
  fRecoWidth[1] = std::abs(minPos[1] - maxPos[1]);
  fRecoWidth[2] = std::abs(minPos[2] - maxPos[2]);

  if(fHitSPE != 0) {
    fRecoPosition[0] /= fHitSPE;
    fRecoPosition[1] /= fHitSPE;
    fRecoPosition[2] /= fHitSPE;
  } else {
    fRecoPosition[0] /= (double)fHitVector.size();
    fRecoPosition[1] /= (double)fHitVector.size();
    fRecoPosition[2] /= (double)fHitVector.size();
  }
  fTimeWidth = std::abs(fLastHitTime-fFirstHitTime);
  fNHit = fHitVector.size();
  fNChannel = channels.size();

  double maxPur=0;
  int gen=-1;
  for (auto& it:fPurity) {
    it.second /= (double)fNHit;
    if (it.second>maxPur) {
      maxPur = it.second;
      gen=it.first;
    }
  }
  //CALL THE CLUSTER MARLEY IF THERE ARE MORE THAN TWO MARLEY HITS IN IT.
  if (gen==1) { fType = 1; }
  else        { fType = 0; }


};

void OpticalCluster::Print() const {
  std::cout << " NChan " << fNChannel;
  std::cout << " NHit " << fNHit;
  std::cout << " Type " << fType;
  std::cout << " SPE " << fHitSPE;
  std::cout << " FirstTime " << fFirstHitTime;
  std::cout << " LastTime " << fLastHitTime;
  std::cout << " TWidth " << fTimeWidth;
  std::cout << " Pos " << fRecoPosition[0] << " " << fRecoPosition[1] << " " << fRecoPosition[2];
  std::cout << " SWidth " << fRecoWidth[0] << " " << fRecoWidth[1] << " " <<fRecoWidth[2] << std::endl;
};
