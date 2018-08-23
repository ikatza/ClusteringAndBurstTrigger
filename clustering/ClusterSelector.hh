#ifndef CLUSTERSELECTOR_HH
#define CLUSTERSELECTOR_HH
#include "WireCluster.hh"


class ClusterSelector {
private:
  int fMinChannel;
  int fMinChannelWidth;
  int fTotalADC;

public:
  ClusterSelector(int MinChannel, int MinChannelWidth, int TotalADC):
    fMinChannel(MinChannel),
    fMinChannelWidth(MinChannelWidth),
    fTotalADC(TotalADC)
    {};
  ClusterSelector():
    fMinChannel(0),
    fMinChannelWidth(0),
    fTotalADC(0)
    {};

  inline bool IsSelected(const WireCluster* clust, int& cut_failed)
    {
      cut_failed=-1;
      bool fail_MinChannel=false, fail_MinChannelWidth=false, fail_TotalADC=false;
      if(clust->GetNChannel() < fMinChannel){
        fail_MinChannel = true;
        cut_failed=0;
      }
    
      if(clust->GetChannelWidth() < fMinChannelWidth){
        fail_MinChannelWidth = true;
        cut_failed=1;
      }

      if(clust->GetHitSumADC() < fTotalADC){
        fail_TotalADC = true;
        cut_failed=2;
      }

      return (fail_MinChannel == false &&
              fail_MinChannelWidth == false &&
              fail_TotalADC == false);
    }

  
  void SetSelectionFlag(WireCluster* clust){
    int flag=-1;
    clust->SetIsSelected(IsSelected(clust, flag));
    clust->SetFailedCut(flag);
  }

  void AreSelected(std::vector<WireCluster*> vec_Clusters)
    {
      for(std::vector<WireCluster*>::iterator it_clust=vec_Clusters.begin();
          it_clust!=vec_Clusters.end(); ++it_clust)
        SetSelectionFlag(*it_clust);
    }

  int  GetMinChannel() const { return fMinChannel;  };
  void SetMinChannel(int inMinChannel=1) { fMinChannel = inMinChannel; };

  int  GetMinChannelWidth() const { return fMinChannelWidth;  };
  void SetMinChannelWidth(int inMinChannelWidth=1) { fMinChannelWidth = inMinChannelWidth; };

  int  GetTotalADC() const { return fTotalADC;  };
  void SetTotalADC(int inTotalADC=1) { fTotalADC = inTotalADC; };

  ~ClusterSelector(){};

public:
//  ClassDef(ClusterSelector,1)  //Simple event class
  
};


#endif
