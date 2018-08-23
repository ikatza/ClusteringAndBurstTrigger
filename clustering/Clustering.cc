#include "Clustering.hh"


//......................................................
int Clustering::ClusterAll(int inNEvent){
  im.SetInputFile(fInputFileName);
  im.SetInputTree(fInputTreeName);
  im.LoadTree();

  if (fERecoXMLFile == "") {
    std::cout << "Not using energy reconstruction!!" << std::endl;
  } else {
    fEReco = new ClusterEnergyEstimator(fERecoXMLFile);
  }
  
  f_Output = new TFile(fOutputFileName.c_str(), "RECREATE");

  t_Output_ClusteredWireHit = new TTree("ClusteredWireHit", "ClusteredWireHit");
  t_Output_ClusteredWireHit->Branch("Cluster",        &out_Cluster,        "Cluster/I"       );
  t_Output_ClusteredWireHit->Branch("Event",          &out_Event,          "Event/I"         );
  t_Output_ClusteredWireHit->Branch("Config",         &out_Config,         "Config/I"        );
  t_Output_ClusteredWireHit->Branch("MarleyIndex",    &out_MarleyIndex,    "MarleyIndex/I"   );
  t_Output_ClusteredWireHit->Branch("StartChan",      &out_StartChan,      "StartChan/I"     );
  t_Output_ClusteredWireHit->Branch("EndChan",        &out_EndChan,        "EndChan/D"       );
  t_Output_ClusteredWireHit->Branch("ChanWidth",      &out_ChanWidth,      "ChanWidth/D"     );
  t_Output_ClusteredWireHit->Branch("NChan",          &out_NChan,          "NChan/D"         );
  t_Output_ClusteredWireHit->Branch("Type",           &out_Type,           "Type/D"          );
  t_Output_ClusteredWireHit->Branch("NHits",          &out_NHits,          "NHits/D"         );
  t_Output_ClusteredWireHit->Branch("SumADC",         &out_SumADC,         "SumADC/D"        );
  t_Output_ClusteredWireHit->Branch("FirstTimeHit",   &out_FirstTimeHit,   "FirstTimeHit/D"  );
  t_Output_ClusteredWireHit->Branch("LastTimeHit",    &out_LastTimeHit,    "LastTimeHit/D"   );
  t_Output_ClusteredWireHit->Branch("TimeWidth",      &out_TimeWidth,      "TimeWidth/D"     );
  t_Output_ClusteredWireHit->Branch("EReco",          &out_EReco,          "EReco/D"         );
  t_Output_ClusteredWireHit->Branch("E_deposited",    &out_E_deposited,    "E_deposited/D"   );
  t_Output_ClusteredWireHit->Branch("NElectron",      &out_NElectron,      "NElectron/D"     );
  t_Output_ClusteredWireHit->Branch("MC_UnderlyingE", &out_MC_UnderlyingE, "MC_UnderlyingE/D");
  t_Output_ClusteredWireHit->Branch("RecClusterPosX", &out_RecClusterPosX, "RecClusterPosX/D");
  t_Output_ClusteredWireHit->Branch("RecClusterPosY", &out_RecClusterPosY, "RecClusterPosY/D");
  t_Output_ClusteredWireHit->Branch("RecClusterPosZ", &out_RecClusterPosZ, "RecClusterPosZ/D");
  t_Output_ClusteredWireHit->Branch("TrClusterPosX",  &out_TrClusterPosX,  "TrClusterPosX/D" );
  t_Output_ClusteredWireHit->Branch("TrClusterPosY",  &out_TrClusterPosY,  "TrClusterPosY/D" );
  t_Output_ClusteredWireHit->Branch("TrClusterPosZ",  &out_TrClusterPosZ,  "TrClusterPosZ/D" );
  t_Output_ClusteredWireHit->Branch("pur_Other"   ,   &out_pur_Other   ,   "pur_Other/D"     );
  t_Output_ClusteredWireHit->Branch("pur_SNnu"    ,   &out_pur_SNnu    ,   "pur_SNnu/D"      );
  t_Output_ClusteredWireHit->Branch("pur_APA"     ,   &out_pur_APA     ,   "pur_APA/D"       );
  t_Output_ClusteredWireHit->Branch("pur_CPA"     ,   &out_pur_CPA     ,   "pur_CPA/D"       );
  t_Output_ClusteredWireHit->Branch("pur_Ar39"    ,   &out_pur_Ar39    ,   "pur_Ar39/D"      );
  t_Output_ClusteredWireHit->Branch("pur_Neutron" ,   &out_pur_Neutron ,   "pur_Neutron/D"   );
  t_Output_ClusteredWireHit->Branch("pur_Krypton" ,   &out_pur_Krypton ,   "pur_Krypton/D"   );
  t_Output_ClusteredWireHit->Branch("pur_Polonium",   &out_pur_Polonium,   "pur_Polonium/D"  );
  t_Output_ClusteredWireHit->Branch("pur_Radon"   ,   &out_pur_Radon   ,   "pur_Radon/D"     );
  t_Output_ClusteredWireHit->Branch("pur_Ar42"    ,   &out_pur_Ar42    ,   "pur_Ar42/D"      );
  t_Output_ClusteredWireHit->Branch("HitView",        &out_HitView);
  t_Output_ClusteredWireHit->Branch("GenType",        &out_GenType);
  t_Output_ClusteredWireHit->Branch("HitChan",        &out_HitChan);
  t_Output_ClusteredWireHit->Branch("HitTime",        &out_HitTime);
  t_Output_ClusteredWireHit->Branch("HitSADC",        &out_HitSADC);
  t_Output_ClusteredWireHit->Branch("HitRMS",         &out_HitRMS );
  t_Output_ClusteredWireHit->Branch("TrPosX",         &out_TrPosX );
  t_Output_ClusteredWireHit->Branch("TrPosY",         &out_TrPosY );
  t_Output_ClusteredWireHit->Branch("TrPosZ",         &out_TrPosZ );

  t_Output_ClusteredOpticalHit = new TTree("ClusteredOpticalHit", "ClusteredOpticalHit");
  t_Output_ClusteredOpticalHit->Branch("Config",         &out_Config,         "Config/I"        );
  t_Output_ClusteredOpticalHit->Branch("Cluster",        &out_Cluster,        "Cluster/I"       );
  t_Output_ClusteredOpticalHit->Branch("MarleyIndex",    &out_MarleyIndex,    "MarleyIndex/I"   );
  t_Output_ClusteredOpticalHit->Branch("Event",          &out_Event,          "Event/I"         );
  t_Output_ClusteredOpticalHit->Branch("YWidth",         &out_YWidth,         "YWidth/D"        );
  t_Output_ClusteredOpticalHit->Branch("ZWidth",         &out_ZWidth,         "ZWidth/D"        );
  t_Output_ClusteredOpticalHit->Branch("NChan",          &out_NChan,          "NChan/D"         );
  t_Output_ClusteredOpticalHit->Branch("Type",           &out_Type,           "Type/D"          );
  t_Output_ClusteredOpticalHit->Branch("NHits",          &out_NHits,          "NHits/D"         );
  t_Output_ClusteredOpticalHit->Branch("SumPE",          &out_SumPE,          "SumPE/D"         );
  t_Output_ClusteredOpticalHit->Branch("FirstTimeHit",   &out_FirstTimeHit,   "FirstTimeHit/D"  );
  t_Output_ClusteredOpticalHit->Branch("LastTimeHit",    &out_LastTimeHit,    "LastTimeHit/D"   );
  t_Output_ClusteredOpticalHit->Branch("TimeWidth",      &out_TimeWidth,      "TimeWidth/D"     );
  t_Output_ClusteredOpticalHit->Branch("RecClusterPosX", &out_RecClusterPosX, "RecClusterPosX/D");
  t_Output_ClusteredOpticalHit->Branch("RecClusterPosY", &out_RecClusterPosY, "RecClusterPosY/D");
  t_Output_ClusteredOpticalHit->Branch("RecClusterPosZ", &out_RecClusterPosZ, "RecClusterPosZ/D");
  t_Output_ClusteredOpticalHit->Branch("Hit_GenType"     , &out_PDSHit_GenType  );
  t_Output_ClusteredOpticalHit->Branch("Hit_X"           , &out_PDSHit_X        );
  t_Output_ClusteredOpticalHit->Branch("Hit_Y"           , &out_PDSHit_Y        );
  t_Output_ClusteredOpticalHit->Branch("Hit_Z"           , &out_PDSHit_Z        );
  t_Output_ClusteredOpticalHit->Branch("Hit_PeakTime"    , &out_PDSHit_PeakTime );
  t_Output_ClusteredOpticalHit->Branch("Hit_Width"       , &out_PDSHit_Width    );
  t_Output_ClusteredOpticalHit->Branch("Hit_PE"          , &out_PDSHit_PE       );
  t_Output_ClusteredOpticalHit->Branch("Hit_OpChannel"   , &out_PDSHit_OpChannel);

  t_Output_TrueInfo = new TTree("TrueInfo", "TrueInfo");
  t_Output_TrueInfo->Branch("Event",    &out_Event, "Event/I");
  t_Output_TrueInfo->Branch("MarlTime", &out_MarlTime);
  t_Output_TrueInfo->Branch("ENu",      &out_ENu     );
  t_Output_TrueInfo->Branch("ENu_Lep",  &out_ENu_Lep );
  t_Output_TrueInfo->Branch("PosX",     &out_PosX    );
  t_Output_TrueInfo->Branch("PosY",     &out_PosY    );
  t_Output_TrueInfo->Branch("PosZ",     &out_PosZ    );
  t_Output_TrueInfo->Branch("PosT",     &out_PosT    );
  t_Output_TrueInfo->Branch("DirX",     &out_DirX    );
  t_Output_TrueInfo->Branch("DirY",     &out_DirY    );
  t_Output_TrueInfo->Branch("DirZ",     &out_DirZ    );

  t_Output_TimingInfo = new TTree("TimingInfo", "TimingInfo");
  t_Output_TimingInfo->Branch("Event",  &out_Event,  "Event/I" );
  t_Output_TimingInfo->Branch("Config", &out_Config, "Config/I");
  t_Output_TimingInfo->Branch("nHit",   &out_NHits,  "nHit/I");
  t_Output_TimingInfo->Branch("TimeOrdering_WireClustTime",         &TimeOrdering_WireClustTime        );
  t_Output_TimingInfo->Branch("SpaceOrdering_WireClustTime",        &SpaceOrdering_WireClustTime       );
  t_Output_TimingInfo->Branch("Clustering_WireClustTime",           &Clustering_WireClustTime          );
  t_Output_TimingInfo->Branch("EnergyReconstruction_WireClustTime", &EnergyReconstruction_WireClustTime);
  t_Output_TimingInfo->Branch("TimeOrdering_OptClustTime",          &TimeOrdering_OptClustTime         );
  t_Output_TimingInfo->Branch("SpaceOrdering_OptClustTime",         &SpaceOrdering_OptClustTime        );
  t_Output_TimingInfo->Branch("Clustering_OptClustTime",            &Clustering_OptClustTime           );
  t_Output_TimingInfo->Branch("EnergyReconstruction_OptClustTime",  &EnergyReconstruction_OptClustTime );
  
  h_ENu_MC      = new TH1D("h_ENu_MC",      "h_ENu_MC",       35,    0, 50  );
  h_MarlTime_MC = new TH1D("h_MarlTime_MC", "h_MarlTime_MC", 100, -0.1, 10.5);
  h_ENu_MC->Sumw2();

  if (inNEvent!=-1) {
    fNEvent = std::min(inNEvent,(int)im.GetEntries());
  } else {
    fNEvent = (int)im.GetEntries();
  }

  TH1I *hNEvents = new TH1I("hNEvents", "hNEvents", 100,0,10e6);
  hNEvents->Fill(fNEvent);
  hNEvents->Write();
  TH1I *hNConfigs = new TH1I("hNConfigs", "hNConfigs", 50,0,50);
  hNConfigs->Fill(fNConfig);
  hNConfigs->Write();

  fClustEng   = new ClusterEngine();
  fClustEng->SetSorting(fSorting);
  fTrigger    = new Trigger();
  fClustSelec = new ClusterSelector();
  fNConfig = (int)fvec_cut_AdjChanTolerance.size();

  std::map<std::pair<int,int>, std::vector<WireHit*>* >     map_unusedHits;
  std::map<std::pair<int,int>, std::vector<WireCluster*>* > map_clusters;
  
  int ConfigBegin = 0;
  int ConfigEnd   = fNConfig;
  if (fConfig != -1) {
    ConfigBegin = fConfig;
    ConfigEnd   = fConfig+1;
  }
  for (size_t iEvent=0; iEvent<fNEvent; ++iEvent) {

    PrintProgress(iEvent,fNEvent);
    im.GetEntry(iEvent);
    out_Event = im.Event;
    bool goodEvent = true;
    for(size_t MarleyEvent=0; MarleyEvent<im.True_MarlSample->size(); ++MarleyEvent) {
      h_ENu_MC     ->Fill(1000*(*im.True_ENu)[MarleyEvent]);
      h_MarlTime_MC->Fill((*im.True_MarlTime)[MarleyEvent]);
      out_MarlTime.push_back((*im.True_MarlTime)[MarleyEvent]);
      out_ENu     .push_back((*im.True_ENu)     [MarleyEvent]);
      out_ENu_Lep .push_back((*im.True_ENu_Lep) [MarleyEvent]);
      out_PosX    .push_back((*im.True_VertX)   [MarleyEvent]);
      out_PosY    .push_back((*im.True_VertY)   [MarleyEvent]);
      out_PosZ    .push_back((*im.True_VertY)   [MarleyEvent]);
      out_PosT    .push_back((*im.True_VertexT) [MarleyEvent]);
      out_DirX    .push_back((*im.True_Dirx)    [MarleyEvent]);
      out_DirY    .push_back((*im.True_Diry)    [MarleyEvent]);
      out_DirZ    .push_back((*im.True_Dirz)    [MarleyEvent]);
    }
    t_Output_TrueInfo->Fill();
    if(!goodEvent) continue;
    //MAKE RECOHIT OBJECTS EVENTWISE FROM THE TREE.
    std::vector<WireHit*>    vec_WireHit;
    std::vector<OpticalHit*> vec_OptHit;

    for(int j = 0; j < im.NColHit; j++) {
      int marley_index=0;
      if (im.Hit_True_MarleyIndex != NULL) {
        marley_index=(*im.Hit_True_MarleyIndex)[j];
      }
      if((*im.Hit_Chan)[j] > fmap_APA_Channel[fNAPA])
        continue;
      
      vec_WireHit.push_back(new WireHit((*im.Hit_View)[j],        (*im.Hit_True_GenType)[j],  (*im.Hit_Chan)[j],
                                        (*im.Hit_Time)[j],        (*im.Hit_SADC)[j],          (*im.Hit_RMS)[j],
                                        (*im.Hit_True_Energy)[j], (*im.Hit_True_EvEnergy)[j], (*im.Hit_True_MainTrID)[j],
                                        0.5*((*im.Hit_X_start)[j]+(*im.Hit_X_end)[j]),
                                        0.5*((*im.Hit_Y_start)[j]+(*im.Hit_Y_end)[j]),
                                        0.5*((*im.Hit_Z_start)[j]+(*im.Hit_Z_end)[j]),
                                        (*im.Hit_True_X)[j],      (*im.Hit_True_Y)[j],        (*im.Hit_True_Z)[j],
                                        marley_index,             (*im.Hit_True_nElec)[j]));
    }
    out_NHits = vec_WireHit.size();
    for(size_t j = 0; j < im.PDS_OpHit_OpChannel->size(); j++) {
      int marley_index=0;
      vec_OptHit.push_back(new OpticalHit((*im.PDS_OpHit_True_GenType)[j],
                                          (*im.PDS_OpHit_X)[j],        (*im.PDS_OpHit_Y)[j],     (*im.PDS_OpHit_Z)[j],
                                          (*im.PDS_OpHit_PeakTime)[j], (*im.PDS_OpHit_Width)[j], (*im.PDS_OpHit_PE)[j],
                                          (*im.PDS_OpHit_OpChannel)[j]));
      vec_OptHit.back()->SetMarleyIndex(marley_index);
    }

    std::random_shuffle(vec_WireHit.begin(), vec_WireHit.end());
    std::random_shuffle(vec_OptHit .begin(), vec_OptHit .end());
    std::random_shuffle(vec_WireHit.begin(), vec_WireHit.end());
    std::random_shuffle(vec_OptHit .begin(), vec_OptHit .end());
    std::random_shuffle(vec_WireHit.begin(), vec_WireHit.end());
    std::random_shuffle(vec_OptHit .begin(), vec_OptHit .end());

    for(fCurrentConfig=ConfigBegin; fCurrentConfig<ConfigEnd; ++fCurrentConfig) {
      fClustEng  ->SetTimeWindow     (fvec_cut_TimeWindowSize  [fCurrentConfig]);
      fClustEng  ->SetTimeWindowOpt  (fvec_cut_TimeWindowOpt[fCurrentConfig]);
      fClustEng  ->SetPositionOpt    (fvec_cut_PositionOpt[fCurrentConfig]);
      fClustEng  ->SetBucketSize     (fvec_cut_BucketSize[fCurrentConfig]);
      fClustEng  ->SetChannelWidth   (fvec_cut_AdjChanTolerance[fCurrentConfig]);
      fClustSelec->SetMinChannel     (fvec_cut_MinChannels     [fCurrentConfig]);
      fClustSelec->SetMinChannelWidth(fvec_cut_MinChanWidth    [fCurrentConfig]);
      fClustSelec->SetTotalADC       (fvec_cut_TotalADC        [fCurrentConfig]);
      fTrigger   ->SetNHitsMin       (fvec_cut_HitsInWindow    [fCurrentConfig]);

      std::vector<WireCluster*>*   vec_Clusters       = new std::vector<WireCluster*>();
      std::vector<WireHit*>*       vec_UnusedHits     = new std::vector<WireHit*>();
      std::vector<OpticalCluster*> vec_OpticalCluster;

      fClustEng->ClusterOpticalHits(vec_OptHit,vec_OpticalCluster);
      FillClusterEngineTimingInfo_Opti(fClustEng);
      FillClusterERecoTimingInfo_Opti (fEReco);

      fClustEng->ClusterHits2(vec_WireHit, vec_Clusters, vec_UnusedHits);
      if (fEReco) fEReco->EstimateEnergy(vec_Clusters);
      FillClusterEngineTimingInfo_Wire(fClustEng);
      FillClusterERecoTimingInfo_Wire (fEReco);

      t_Output_TimingInfo->Fill();

      for(unsigned int k = 0; k < vec_Clusters->size(); k++) {
        WireCluster* aCluster = vec_Clusters->at(k);
        fClustSelec->SetSelectionFlag  (aCluster);
        fTrigger   ->SetTriggerStateFor(aCluster);
      }

      for(size_t i=0; i<vec_UnusedHits->size();++i)
        FillUnusedHitInfo((*vec_UnusedHits)[i]);

      for(size_t i=0; i<vec_Clusters->size();++i) {
        FillClusterInfo((*vec_Clusters)[i]);
        delete (*vec_Clusters)[i];
        (*vec_Clusters)[i]=NULL;
      }

      if (fPrintLevel > -1) {
        std::cout << "----------------------------------------------" << std::endl;
        std::cout << "Config           " << fCurrentConfig            << std::endl;
        std::cout << "nWireClusters    " << vec_Clusters->size()      << std::endl;
        std::cout << "nOpticalClusters " << vec_OpticalCluster.size() << std::endl;
      }

      std::sort(vec_OpticalCluster.begin(), vec_OpticalCluster.end(), OpticalClusterOrderedInTimePtr);

      for(size_t i=0; i<vec_OpticalCluster.size();++i) {
        FillClusterInfo(vec_OpticalCluster[i]);
        delete vec_OpticalCluster[i];
        vec_OpticalCluster[i]=NULL;
      }
    }
      
    for(size_t j = 0; j < vec_WireHit.size(); j++) {
      delete vec_WireHit[j];
      vec_WireHit[j] = NULL;
    }
    for(size_t j = 0; j < vec_OptHit.size(); j++) {
      delete vec_OptHit[j];
      vec_OptHit[j] = NULL;
    }
    vec_WireHit.clear();
  }

  std::cout << "Writing all the trees and closing." << std::endl;
  h_ENu_MC     ->Write();
  h_MarlTime_MC->Write();

  t_Output_ClusteredWireHit   ->Write();
  t_Output_ClusteredOpticalHit->Write();
  t_Output_TrueInfo           ->Write();
  t_Output_TimingInfo         ->Write();
  
  return 0;

};

void Clustering::FillClusterEngineTimingInfo_Opti(ClusterEngine* clusterEngine) {
  TimeOrdering_OptClustTime .clear();
  SpaceOrdering_OptClustTime.clear();
  Clustering_OptClustTime   .clear();
  out_Config = fCurrentConfig;
  out_Event  = im.Event;
  TimeOrdering_OptClustTime  = clusterEngine->GetElapsedTime_TimeOrdering ();
  SpaceOrdering_OptClustTime = clusterEngine->GetElapsedTime_SpaceOrdering();
  Clustering_OptClustTime    = clusterEngine->GetElapsedTime_Clustering   ();
  
};

void Clustering::FillClusterEngineTimingInfo_Wire(ClusterEngine* clusterEngine) {
  TimeOrdering_WireClustTime .clear();
  SpaceOrdering_WireClustTime.clear();
  Clustering_WireClustTime   .clear();
  out_Config = fCurrentConfig;
  out_Event  = im.Event;
  TimeOrdering_WireClustTime  = clusterEngine->GetElapsedTime_TimeOrdering ();
  SpaceOrdering_WireClustTime = clusterEngine->GetElapsedTime_SpaceOrdering();
  Clustering_WireClustTime    = clusterEngine->GetElapsedTime_Clustering   ();
  
};

void Clustering::FillClusterERecoTimingInfo_Wire(ClusterEnergyEstimator* clusterEReco) {
  EnergyReconstruction_WireClustTime.clear();
  EnergyReconstruction_WireClustTime = clusterEReco->GetElapsedTime();
};

void Clustering::FillClusterERecoTimingInfo_Opti(ClusterEnergyEstimator* clusterEReco) {
  EnergyReconstruction_OptClustTime = {0};
};

void Clustering::FillClusterInfo(OpticalCluster* clust) {
  ResetFillVariable();
  //FILL THE OUTPUT TREE.
  out_Config         = fCurrentConfig;
  out_Cluster        = fvec_OptClusterCount[0];
  out_Event          = im.Event;
  out_MarleyIndex    = clust->GetMarleyIndex();
  out_YWidth         = clust->GetRecoWidth(1);
  out_ZWidth         = clust->GetRecoWidth(2);
  out_NChan          = clust->GetNChannel();
  out_Type           = clust->GetType();
  out_NHits          = clust->GetNHit();
  out_SumPE          = clust->GetHitSPE();
  out_FirstTimeHit   = clust->GetFirstHitTime();
  out_LastTimeHit    = clust->GetLastHitTime();
  out_TimeWidth      = clust->GetTimeWidth();
  out_RecClusterPosX = clust->GetRecoPosition(0);
  out_RecClusterPosY = clust->GetRecoPosition(1);
  out_RecClusterPosZ = clust->GetRecoPosition(2);

  if(fPrintLevel == 1 || fPrintLevel > 2)
    clust->Print();

  for(auto const& hit : clust->GetHitVector()) {
    out_PDSHit_GenType  .push_back(hit->GetGenType());
    out_PDSHit_X        .push_back(hit->GetRecoPosition(0));
    out_PDSHit_Y        .push_back(hit->GetRecoPosition(1));
    out_PDSHit_Z        .push_back(hit->GetRecoPosition(2));
    out_PDSHit_PeakTime .push_back(hit->GetTime   ());
    out_PDSHit_Width    .push_back(hit->GetWidth  ());
    out_PDSHit_PE       .push_back(hit->GetSPE    ());
    out_PDSHit_OpChannel.push_back(hit->GetChannel());
  }
  t_Output_ClusteredOpticalHit->Fill();

};

void Clustering::FillClusterInfo(WireCluster* clust) {
  
  ResetFillVariable();
  //FILL THE OUTPUT TREE.
  out_Event = im.Event;
  out_Config         = fCurrentConfig;
  out_Cluster        = fvec_ClusterCount[fCurrentConfig];
  out_MarleyIndex    = clust->GetMarleyIndex();
  out_StartChan      = clust->GetStartChannel();
  out_EndChan        = clust->GetEndChannel();
  out_ChanWidth      = clust->GetChannelWidth();
  out_NChan          = clust->GetNChannel();
  out_Type           = clust->GetType();
  out_NHits          = clust->GetNHit();
  out_SumADC         = clust->GetHitSumADC();
  out_FirstTimeHit   = clust->GetFirstTimeHit();
  out_LastTimeHit    = clust->GetLastTimeHit();
  out_TimeWidth      = clust->GetTimeWidth();
  out_EReco          = clust->GetEReco();
  out_MC_UnderlyingE = clust->GetMC_UnderlyingE();
  out_NElectron      = clust->GetNElectron();
  out_TrClusterPosX  = clust->GetTruePosition(0);
  out_TrClusterPosY  = clust->GetTruePosition(1);
  out_TrClusterPosZ  = clust->GetTruePosition(2);
  out_RecClusterPosX = clust->GetRecoPosition(0);
  out_RecClusterPosY = clust->GetRecoPosition(1);
  out_RecClusterPosZ = clust->GetRecoPosition(2);
  out_pur_Other      = clust->GetPurity(kOther   );
  out_pur_SNnu       = clust->GetPurity(kSNnu    );
  out_pur_APA        = clust->GetPurity(kAPA     );
  out_pur_CPA        = clust->GetPurity(kCPA     );
  out_pur_Ar39       = clust->GetPurity(kAr39    );
  out_pur_Neutron    = clust->GetPurity(kNeutron );
  out_pur_Krypton    = clust->GetPurity(kKrypton );
  out_pur_Polonium   = clust->GetPurity(kPolonium);
  out_pur_Radon      = clust->GetPurity(kRadon   );
  out_pur_Ar42       = clust->GetPurity(kAr42    );
  if(fPrintLevel > 1)
    clust->Print();

  std::vector<WireHit*> hits = clust->GetHitVector();
  out_E_deposited=0.;
  for(auto const& hit : clust->GetHitVector())
  {
    out_HitView.push_back(hit->GetHitView());
    out_GenType.push_back(hit->GetGenType());
    out_HitChan.push_back(hit->GetHitChan());
    out_HitTime.push_back(hit->GetHitTime());
    out_HitSADC.push_back(hit->GetHitSADC());
    out_HitRMS .push_back(hit->GetHitRMS ());
    out_TrPosX .push_back(hit->GetTruePosition(0));
    out_TrPosY .push_back(hit->GetTruePosition(1));
    out_TrPosZ .push_back(hit->GetTruePosition(2));
    out_E_deposited += hit->GetTrDepositedE();
  }
  if(clust->GetTriggerFlag() && clust->GetIsSelected()){
    // std::cout << "EVENT " << out_Event << std::endl;
    // clust->Print();
    t_Output_ClusteredWireHit->Fill();
    fvec_ClusterCount[fCurrentConfig]++;
  }

};

void Clustering::ResetFillVariable(){
  out_Cluster        = 0;
  out_Event          = 0;
  out_Config         = 0;
  out_StartChan      = 0;
  out_EndChan        = 0;
  out_ChanWidth      = 0;
  out_NChan          = 0;
  out_Type           = 0;
  out_NHits          = 0;
  out_SumADC         = 0;
  out_FirstTimeHit   = 0;
  out_LastTimeHit    = 0;
  out_TimeWidth      = 0;
  out_E_deposited    = 0;
  out_NElectron      = 0;
  out_MC_UnderlyingE = 0;
  out_TrClusterPosX  = 0;
  out_TrClusterPosY  = 0;
  out_TrClusterPosZ  = 0;
  out_RecClusterPosX = 0;
  out_RecClusterPosY = 0;
  out_RecClusterPosZ = 0;
  out_pur_Other      = 0;
  out_pur_SNnu       = 0;
  out_pur_APA        = 0;
  out_pur_CPA        = 0;
  out_pur_Ar39       = 0;
  out_pur_Neutron    = 0;
  out_pur_Krypton    = 0;
  out_pur_Polonium   = 0;
  out_pur_Radon      = 0;
  out_pur_Ar42       = 0;
  out_YWidth         = 0;
  out_ZWidth         = 0;
  out_SumPE          = 0;
  out_EReco          = 0;
  out_MarleyIndex    =-1;
  out_MarlTime.clear();
  out_ENu     .clear();
  out_ENu_Lep .clear();
  out_PosX    .clear();
  out_PosY    .clear();
  out_PosZ    .clear();
  out_PosT    .clear();
  out_DirX    .clear();
  out_DirY    .clear();
  out_DirZ    .clear();

  out_HitView .clear();
  out_GenType .clear();
  out_HitChan .clear();
  out_HitTime .clear();
  out_HitSADC .clear();
  out_HitSPE  .clear();
  out_HitRMS  .clear();
  out_HitWidth.clear();
  out_TrPosX  .clear();
  out_TrPosY  .clear();
  out_TrPosZ  .clear();
  out_PDSHit_GenType  .clear();
  out_PDSHit_X        .clear();
  out_PDSHit_Y        .clear();
  out_PDSHit_Z        .clear();
  out_PDSHit_PeakTime .clear();
  out_PDSHit_Width    .clear();
  out_PDSHit_PE       .clear();
  out_PDSHit_OpChannel.clear();

  TimeOrdering_WireClustTime        .clear();
  SpaceOrdering_WireClustTime       .clear();
  Clustering_WireClustTime          .clear();
  EnergyReconstruction_WireClustTime.clear();
  TimeOrdering_OptClustTime         .clear();
  SpaceOrdering_OptClustTime        .clear();
  Clustering_OptClustTime           .clear();
  EnergyReconstruction_OptClustTime .clear();

};
