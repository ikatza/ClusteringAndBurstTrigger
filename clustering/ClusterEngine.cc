#include "ClusterEngine.hh"

//ClassImp(ClusterEngine)

void ClusterEngine::ClusterHits2_stable_sort(std::vector<WireHit*>& vec_Hits,
                                 std::vector<WireCluster*>* vec_Clusters,
                                 std::vector<WireHit*>* vec_UnusedHits) {

  ResetTime();
  fStopwatch0 = clock();
  std::stable_sort(vec_Hits.begin(), vec_Hits.end(), WireHitOrderedInSpacePtr);
  SpaceOrdering_Time.push_back(((float)(clock()-fStopwatch0))/((float)CLOCKS_PER_SEC));
   
  vec_Clusters->clear();
  (void) vec_UnusedHits;

  if(vec_Hits.size() < 2) return;
  for(size_t i = 0; i < vec_Hits.size()-1; i++){

    std::vector<WireHit*> vec_TempHits;
    if(std::abs(vec_Hits[i]->GetHitChan()-vec_Hits[i+1]->GetHitChan())<=fChannelWidth)
    {
      int channelCount = 1;

      vec_TempHits.push_back(vec_Hits[i  ]);
      vec_TempHits.push_back(vec_Hits[i+1]);

      while((i+channelCount+1)<vec_Hits.size() &&
            std::abs(vec_Hits[i+channelCount]->GetHitChan()-vec_Hits[i+channelCount+1]->GetHitChan()) <= fChannelWidth)
      {
        vec_TempHits.push_back(vec_Hits[i + channelCount + 1]);
        channelCount++;
      }

      i = i + channelCount;
      fStopwatch1 = clock();
      std::stable_sort(vec_TempHits.begin(), vec_TempHits.end(), WireHitOrderedInTimePtr);
      TimeOrdering_Time.push_back((float(clock() - fStopwatch1))/((float)CLOCKS_PER_SEC));

      for(size_t j = 0; j < vec_TempHits.size()-1; j++)
      {
        std::vector<WireHit*> vec_TempHits2;
        if(std::abs(vec_TempHits[j]->GetHitTime()-vec_TempHits[j+1]->GetHitTime())<=fTimeWindow)
        {
          int timeCount = 1;

          vec_TempHits2.push_back(vec_TempHits[j  ]);
          vec_TempHits2.push_back(vec_TempHits[j+1]);

          while((j+timeCount+1)<vec_TempHits.size() &&
                std::abs(vec_TempHits[j+timeCount]->GetHitTime()-vec_TempHits[j+timeCount+1]->GetHitTime())
                <=fTimeWindow)
          {
            if(vec_TempHits[j + timeCount + 1]->GetHitSADC()>fMinHitADC)
              vec_TempHits2.push_back(vec_TempHits[j + timeCount + 1]);
            timeCount++;
          }

          j = j + timeCount;
          vec_Clusters->push_back(new WireCluster(vec_TempHits2));
          vec_TempHits2.clear();
          //fVecClusters.push_back(temp);
        }
      }
      vec_TempHits.clear();
    }
  }
  Clustering_Time.push_back(((float)(clock() - fStopwatch0))/((float)CLOCKS_PER_SEC));
  
};


void ClusterEngine::ClusterHits2_sort(std::vector<WireHit*>& vec_Hits,
                                 std::vector<WireCluster*>* vec_Clusters,
                                 std::vector<WireHit*>* vec_UnusedHits) {

  ResetTime();
  fStopwatch0 = clock();
  std::sort(vec_Hits.begin(), vec_Hits.end(), WireHitOrderedInSpacePtr);
  SpaceOrdering_Time.push_back(((float)(clock()-fStopwatch0))/((float)CLOCKS_PER_SEC));
   
  vec_Clusters->clear();
  (void) vec_UnusedHits;
  
  if(vec_Hits.size() < 2) return;
  for(size_t i = 0; i < vec_Hits.size()-1; i++){
    
    std::vector<WireHit*> vec_TempHits;
    if(std::abs(vec_Hits[i]->GetHitChan()-vec_Hits[i+1]->GetHitChan())<=fChannelWidth)
    {
      int channelCount = 1;
      
      vec_TempHits.push_back(vec_Hits[i  ]);
      vec_TempHits.push_back(vec_Hits[i+1]);

      while((i+channelCount+1)<vec_Hits.size() &&
            std::abs(vec_Hits[i+channelCount]->GetHitChan()-vec_Hits[i+channelCount+1]->GetHitChan()) <= fChannelWidth)
      {
        vec_TempHits.push_back(vec_Hits[i + channelCount + 1]);
        channelCount++;
      }

      i = i + channelCount;
      fStopwatch1 = clock();
      std::sort(vec_TempHits.begin(), vec_TempHits.end(), WireHitOrderedInTimePtr);
      TimeOrdering_Time.push_back((float(clock() - fStopwatch1))/((float)CLOCKS_PER_SEC));

      for(size_t j = 0; j < vec_TempHits.size()-1; j++)
      {
        std::vector<WireHit*> vec_TempHits2;
        if(std::abs(vec_TempHits[j]->GetHitTime()-vec_TempHits[j+1]->GetHitTime())<=fTimeWindow)
        {
          int timeCount = 1;
                                     
          vec_TempHits2.push_back(vec_TempHits[j  ]);
          vec_TempHits2.push_back(vec_TempHits[j+1]);

          while((j+timeCount+1)<vec_TempHits.size() && 
                std::abs(vec_TempHits[j+timeCount]->GetHitTime()-vec_TempHits[j+timeCount+1]->GetHitTime())
                <=fTimeWindow)
          {
            if(vec_TempHits[j + timeCount + 1]->GetHitSADC()>fMinHitADC)
              vec_TempHits2.push_back(vec_TempHits[j + timeCount + 1]);
            timeCount++;
          }

          j = j + timeCount;
          vec_Clusters->push_back(new WireCluster(vec_TempHits2));
          vec_TempHits2.clear();
          //fVecClusters.push_back(temp);
        }
      }
      vec_TempHits.clear();
    }
  }
  Clustering_Time.push_back(((float)(clock() - fStopwatch0))/((float)CLOCKS_PER_SEC));
  
};


void ClusterEngine::ClusterOpticalHits(std::vector<OpticalHit*>& vec_OptHit,
                                       std::vector<OpticalCluster*>& vec_OptCluster) {

  ResetTime();
  fStopwatch0 = clock();
  std::sort(vec_OptHit.begin(), vec_OptHit.end(), OpticalHitOrderedInTimePtr);
  TimeOrdering_Time.push_back(((float)(clock() - fStopwatch0))/((float)CLOCKS_PER_SEC));
  vec_OptCluster.clear();

  if(vec_OptHit.size() < 2) return;
  for(size_t i = 0; i < vec_OptHit.size()-1; i++){

    std::vector<OpticalHit*> vec_TempHit;
    if(std::abs(vec_OptHit[i]->GetTime()-vec_OptHit[i+1]->GetTime())<=fTimeWindowOpt)
    {
      int channelCount = 1;

      vec_TempHit.push_back(vec_OptHit[i  ]);
      vec_TempHit.push_back(vec_OptHit[i+1]);

      while((i+channelCount+1)<vec_OptHit.size() &&
            std::abs(vec_OptHit[i+channelCount]->GetTime()-vec_OptHit[i+channelCount+1]->GetTime())<=fTimeWindowOpt)
      {
        vec_TempHit.push_back(vec_OptHit[i + channelCount + 1]);
        channelCount++;
      }

      i = i + channelCount;
      std::map<OpticalHit*, std::vector<OpticalHit*>> bucket;
      for (size_t ii = 0; ii < vec_TempHit.size(); ii++){
          // bucket[vec_TempHit[ii]].push_back(vec_TempHit[0]);
          for(size_t j = ii; j < vec_TempHit.size(); j++){
              double time_diff = vec_TempHit[j]->GetTime() -  vec_TempHit[ii]->GetTime();
              if (time_diff < fBucketSize)
                  bucket[vec_TempHit[ii]].push_back(vec_TempHit[j]);
          }
      }
      int max_threshold = 0;
      OpticalHit* optical_hit_idx;
      for (size_t iii=0; iii<vec_TempHit.size(); iii++){
          size_t bucket_size = bucket[vec_TempHit[iii]].size();
          if (bucket_size > max_threshold){
              max_threshold = bucket_size;
              optical_hit_idx = vec_TempHit[iii];
          }
      }

      //std::sort(vec_TempHit.begin(), vec_TempHit.end(), OpticalHitOrderedInSpacePtr);
      //if(bucket[optical_hit_idx].size() > 4 )
      // std::cout << "vec_TempHit " << vec_TempHit.size() << "\t" ;
      // std::cout << "bucket size " << bucket[optical_hit_idx].size() << std::endl;

      fStopwatch1 = clock();
      std::sort(vec_TempHit.begin(), vec_TempHit.end(), OpticalHitOrderedInSpacePtr);
      SpaceOrdering_Time.push_back(((float)(clock() - fStopwatch1))/((float)CLOCKS_PER_SEC));
      for(size_t j = 0; j < vec_TempHit.size()-1; j++){
        
        std::vector<OpticalHit*> vec_TempHit2;
        if(std::abs(bucket[optical_hit_idx][j]->GetRecoPosition(2)-bucket[optical_hit_idx][j+1]->GetRecoPosition(2))<=fPositionOpt)
        {
          int PosCount = 1;

          vec_TempHit2.push_back(bucket[optical_hit_idx][j  ]);
          vec_TempHit2.push_back(bucket[optical_hit_idx][j+1]);

          while((j+PosCount+1)<bucket[optical_hit_idx].size() &&
                std::abs(bucket[optical_hit_idx][j+PosCount]->GetRecoPosition(2)-bucket[optical_hit_idx][j+PosCount+1]->GetRecoPosition(2))<=fPositionOpt)
          {
            vec_TempHit2.push_back(bucket[optical_hit_idx][j + PosCount + 1]);
            PosCount++;
          }

          j = j + PosCount;
          vec_OptCluster.push_back(new OpticalCluster(vec_TempHit2));
        }
      }
    }
  }
  Clustering_Time.push_back(((float)(clock() - fStopwatch0))/((float)CLOCKS_PER_SEC));

};


void ClusterEngine::ClusterOpticalHits2(std::vector<OpticalHit*>& vec_OptHit,
                                       std::vector<OpticalCluster*>& vec_OptCluster) {

  std::sort(vec_OptHit.begin(), vec_OptHit.end(), OpticalHitOrderedInTimePtr);
  vec_OptCluster.clear();

  if(vec_OptHit.size() < 2) return;
  for(size_t i = 0; i < vec_OptHit.size()-1; i++){

    std::vector<OpticalHit*> vec_TempHit;
    if(std::abs(vec_OptHit[i]->GetTime()-vec_OptHit[i+1]->GetTime())<=fTimeWindowOpt)
    {
      int channelCount = 1;

      vec_TempHit.push_back(vec_OptHit[i  ]);
      vec_TempHit.push_back(vec_OptHit[i+1]);

      while((i+channelCount+1)<vec_OptHit.size() &&
            std::abs(vec_OptHit[i+channelCount]->GetTime()-vec_OptHit[i+channelCount+1]->GetTime())<=fTimeWindowOpt)
      {
        vec_TempHit.push_back(vec_OptHit[i + channelCount + 1]);
        channelCount++;
      }

      i = i + channelCount;
      std::sort(vec_TempHit.begin(), vec_TempHit.end(), OpticalHitOrderedInSpacePtr);
      for(size_t j = 0; j < vec_TempHit.size()-1; j++){

        std::vector<OpticalHit*> vec_TempHit2;
        if(std::abs(vec_TempHit[j]->GetRecoPosition(2)-vec_TempHit[j+1]->GetRecoPosition(2))<=fPositionOpt)
        {
          int PosCount = 1;

          vec_TempHit2.push_back(vec_TempHit[j  ]);
          vec_TempHit2.push_back(vec_TempHit[j+1]);

          while((j+PosCount+1)<vec_TempHit.size() &&
                std::abs(vec_TempHit[j+PosCount]->GetRecoPosition(2)-vec_TempHit[j+PosCount+1]->GetRecoPosition(2))<=fPositionOpt)
          {
            vec_TempHit2.push_back(vec_TempHit[j + PosCount + 1]);
            PosCount++;
          }

          j = j + PosCount;
          vec_OptCluster.push_back(new OpticalCluster(vec_TempHit2));
        }
      }
    }
  }
};


void ClusterEngine::ClusterHits3(const std::vector<WireHit*>& vec_Hits,
                                 std::vector<WireCluster*>* vec_Clusters,
                                 std::vector<WireHit*>* vec_UnusedHits)
{
  vec_Clusters->clear();
  (void) vec_UnusedHits;
  if(vec_Hits.size() < 2) return;
  for(size_t i = 0; i < vec_Hits.size()-1; i++){

    std::vector<WireHit*> vec_TempHits;
    if(std::abs(vec_Hits[i]->GetHitChan()-vec_Hits[i+1]->GetHitChan())<=fChannelWidth &&
       std::abs(vec_Hits[i]->GetHitTime()-vec_Hits[i+1]->GetHitTime())<=fTimeWindow)
    {
      int channelCount = 1;

      vec_TempHits.push_back(vec_Hits[i  ]);
      vec_TempHits.push_back(vec_Hits[i+1]);

      while((i+channelCount+1)<vec_Hits.size() &&
            std::abs(vec_Hits[i+channelCount]->GetHitChan()-vec_Hits[i+channelCount+1]->GetHitChan()) <=fChannelWidth &&
            std::abs(vec_Hits[i+channelCount]->GetHitTime()-vec_Hits[i+channelCount+1]->GetHitTime()) <= fTimeWindow)
      {
        vec_TempHits.push_back(vec_Hits[i + channelCount + 1]);
        channelCount++;
      }

      i = i + channelCount;
      vec_Clusters->push_back(new WireCluster(vec_TempHits));
      vec_TempHits.clear();
    }
  }

};



void ClusterEngine::ClusterHits(const std::vector<WireHit*>& vec_Hits,
                                std::vector<WireCluster*>* vec_Clusters,
                                std::vector<WireHit*>* vec_UnusedHits)
{
  vec_Clusters->clear();
  vec_UnusedHits->clear();

  std::vector<WireHit*> LocalHits = vec_Hits;
  std::sort(LocalHits.begin(), LocalHits.end(), WireHitOrderedInSpacePtr);

  for(std::vector<WireHit*>::iterator it_hit1 = LocalHits.begin();
      it_hit1 != LocalHits.end()-1;++it_hit1)
  {
    std::vector<WireHit*> vec_ClusteredHits;
    vec_ClusteredHits.clear();

    WireHit* hit1 = *it_hit1;
    bool hit1InACluster = false;
    int i = 0;
    for(std::vector<WireHit*>::iterator it_hit2 = it_hit1+1;
        it_hit2 != LocalHits.end();++it_hit2)
    {
      //
      i++;
      WireHit* hit2 = *it_hit2;
      double chanDiff = std::abs(hit1->GetHitChan()-
                                 hit2->GetHitChan());
      double timeDiff = std::abs(hit1->GetHitTime()-
                                 hit2->GetHitTime());

      if(chanDiff<fChannelWidth &&
         timeDiff<fTimeWindow   &&
         std::find(vec_ClusteredHits.begin(), vec_ClusteredHits.end(), hit2) == vec_ClusteredHits.end())
      {
        hit1InACluster = true;

        vec_ClusteredHits.push_back(hit2);
      }
    }

    if(hit1InACluster &&
       std::find(vec_ClusteredHits.begin(), vec_ClusteredHits.end(), hit1) == vec_ClusteredHits.end()){
      vec_ClusteredHits.push_back(hit1);
      vec_Clusters->push_back(new WireCluster(vec_ClusteredHits));
    }else{
      vec_UnusedHits->push_back(hit1);
    }

  }

};
