#ifndef TRIGGER_HH
#define TRIGGER_HH
#include "WireCluster.hh"

class Trigger {
public:
  Trigger(int inNHitsMin): fNHitsMin(inNHitsMin){};
  Trigger(): fNHitsMin(0){};

  bool GetTriggerStateFor(WireCluster* clust){ return goodClusterHits(clust, fNHitsMin); }
  void SetTriggerStateFor(WireCluster* clust)
    {
      clust->SetTriggerFlag(goodClusterHits(clust, fNHitsMin));
    };
  void SetTriggerStateFor(std::vector<WireCluster*> clusters)
    {
      for(std::vector<WireCluster*>::iterator it_clust=clusters.begin();
          it_clust!=clusters.end(); ++it_clust)
        (*it_clust)->SetTriggerFlag(goodClusterHits((*it_clust), fNHitsMin));
    };

  int  GetNHitsMin() const { return fNHitsMin;  };
  void SetNHitsMin(int inNHitsMin=1) { fNHitsMin = inNHitsMin; };
  ~Trigger(){};

private:
  int fNHitsMin;

public:
//  ClassDef(Trigger,1)  //Simple event class

};


#endif
