#ifndef INPUTMANAGER_H
#define INPUTMANAGER_H

#include "TFile.h"
#include "TTree.h"

#include <vector>
#include <iostream>

class InputManager{
private:
  TFile* f_Input;
  std::string treename;
  std::string filename;
public:
  TTree* t_Input;

  int Run;
  int SubRun;
  int Event;
  
  int NTotHit    ;
  int NColHit    ;
  int NIndHit    ;
  int NHitNoBT   ;
  int NFlash     ;
  int NFlashNoBT ;
  int TotGen_Marl;
  int TotGen_APA ;
  int TotGen_CPA ;
  int TotGen_Ar39;
  int TotGen_Neut;
  int TotGen_Kryp;
  int TotGen_Plon;
  int TotGen_Rdon;
  int TotGen_Ar42;   
 
  std::vector<int>                 * Hit_View                 ;
  std::vector<int>                 * Hit_Size                 ;
  std::vector<int>                 * Hit_TPC                  ;
  std::vector<int>                 * Hit_Chan                 ;
  std::vector<int>                 * Hit_X_start              ;
  std::vector<int>                 * Hit_Y_start              ;
  std::vector<int>                 * Hit_Z_start              ;
  std::vector<int>                 * Hit_X_end                ;
  std::vector<int>                 * Hit_Y_end                ;
  std::vector<int>                 * Hit_Z_end                ;
  std::vector<float>               * Hit_Time                 ;
  std::vector<float>               * Hit_RMS                  ;
  std::vector<float>               * Hit_SADC                 ;
  std::vector<float>               * Hit_Int                  ;
  std::vector<float>               * Hit_Peak                 ;
  std::vector<int>                 * Hit_True_GenType         ;
  std::vector<int>                 * Hit_True_MainTrID        ;
  std::vector<float>               * Hit_True_EvEnergy        ;
  std::vector<float>               * Hit_True_X               ;                  
  std::vector<float>               * Hit_True_Y               ;                  
  std::vector<float>               * Hit_True_Z               ;                  
  std::vector<float>               * Hit_True_Energy          ;             
  std::vector<float>               * Hit_True_nElec           ;
  std::vector<int>                 * Hit_True_nIDEs           ; 
  std::vector<int>                 * Hit_True_MarleyIndex     ;
  
  std::vector<int>                 * Hit_AdjM5SADC            ;
  std::vector<int>                 * Hit_AdjM2SADC            ;
  std::vector<int>                 * Hit_AdjM1SADC            ;
  std::vector<int>                 * Hit_AdjP1SADC            ;
  std::vector<int>                 * Hit_AdjP2SADC            ;
  std::vector<int>                 * Hit_AdjP5SADC            ;
  std::vector<int>                 * Hit_AdjM5Chan            ;
  std::vector<int>                 * Hit_AdjM2Chan            ;
  std::vector<int>                 * Hit_AdjM1Chan            ;
  std::vector<int>                 * Hit_AdjP1Chan            ;
  std::vector<int>                 * Hit_AdjP2Chan            ;
  std::vector<int>                 * Hit_AdjP5Chan            ;

  std::vector<int>                 * PDS_Flash_FlashID        ;
  std::vector<float>               * PDS_Flash_YCenter        ;
  std::vector<float>               * PDS_Flash_ZCenter        ;
  std::vector<float>               * PDS_Flash_YWidth         ;
  std::vector<float>               * PDS_Flash_ZWidth         ;
  std::vector<float>               * PDS_Flash_Time           ;
  std::vector<float>               * PDS_Flash_TimeWidth      ;
  std::vector<float>               * PDS_Flash_TotalPE        ;
  std::vector<float>               * PDS_Flash_True_Distance  ;
  std::vector<float>               * PDS_Flash_True_GenType   ;

  std::vector<int> 	           * PDS_OpHit_OpChannel      ;
  std::vector<double>              * PDS_OpHit_X              ;
  std::vector<double> 	           * PDS_OpHit_Y              ;
  std::vector<double> 	           * PDS_OpHit_Z              ;
  std::vector<double> 	           * PDS_OpHit_PeakTimeAbs    ;
  std::vector<double> 	           * PDS_OpHit_PeakTime       ;
  std::vector<unsigned short>      * PDS_OpHit_Frame          ;
  std::vector<double>              * PDS_OpHit_Width          ;
  std::vector<double> 	           * PDS_OpHit_Area           ;
  std::vector<double> 	           * PDS_OpHit_Amplitude      ;
  std::vector<double> 	           * PDS_OpHit_PE             ;
  std::vector<double> 	           * PDS_OpHit_FastToTotal    ;
  std::vector<int>                 * PDS_OpHit_True_GenType   ;
  std::vector<double>              * PDS_OpHit_True_Energy    ;
  std::vector<int>                 * PDS_OpHit_True_TrackID   ;
  std::vector<int>                 * PDS_OpHit_True_GenTypeAll;
  std::vector<double>              * PDS_OpHit_True_EnergyAll ;
  std::vector<int>                 * PDS_OpHit_True_TrackIDAll;
  std::vector<int>                 * PDS_OpHit_True_IndexAll  ;

  std::vector<int>                 * True_VertexChan          ;
  std::vector<int>                 * True_Nu_Type             ;
  std::vector<int>                 * True_Nu_Lep_Type         ;
  std::vector<int>                 * True_Mode                ;
  std::vector<int>                 * True_CCNC                ;
  std::vector<int>                 * True_HitNucleon          ;
  std::vector<int>                 * True_Target              ;
  std::vector<int>                 * True_MarlSample          ;
  std::vector<float>               * True_MarlTime            ;
  std::vector<float>               * True_MarlWeight          ;
  std::vector<float>               * True_ENu                 ;
  std::vector<float>               * True_ENu_Lep             ;
  std::vector<float>               * True_VertX               ;
  std::vector<float>               * True_VertY               ;
  std::vector<float>               * True_VertZ               ;
  std::vector<float>               * True_VertexT             ;
  std::vector<float>               * True_Px                  ;
  std::vector<float>               * True_Py                  ;
  std::vector<float>               * True_Pz                  ;
  std::vector<float>               * True_Dirx                ;
  std::vector<float>               * True_Diry                ;
  std::vector<float>               * True_Dirz                ;
  std::vector<float>               * True_Time                ;
                                                              
  std::vector<int>                 * True_Bck_Mode            ;
  std::vector<double>              * True_Bck_VertX           ;
  std::vector<double>              * True_Bck_VertY           ;
  std::vector<double>              * True_Bck_VertZ           ;
  std::vector<double>              * True_Bck_Time            ;
  std::vector<double>              * True_Bck_Energy          ;
  std::vector<int>                 * True_Bck_PDG             ;

  std::vector<int>                 * True_Bck_ModeAll         ;
  std::vector<double>              * True_Bck_VertXAll        ;
  std::vector<double>              * True_Bck_VertYAll        ;
  std::vector<double>              * True_Bck_VertZAll        ;
  std::vector<double>              * True_Bck_TimeAll         ;
  std::vector<double>              * True_Bck_EndVertXAll     ;
  std::vector<double>              * True_Bck_EndVertYAll     ;
  std::vector<double>              * True_Bck_EndVertZAll     ;
  std::vector<double>              * True_Bck_EndTimeAll      ;
  std::vector<double>              * True_Bck_EnergyAll       ;
  std::vector<int>                 * True_Bck_PDGAll          ;
  std::vector<int>                 * True_Bck_IDAll           ;

  
  ~InputManager(){
    delete Hit_View                 ; Hit_View                  = NULL;
    delete Hit_Size                 ; Hit_Size                  = NULL;
    delete Hit_TPC                  ; Hit_TPC                   = NULL;
    delete Hit_Chan                 ; Hit_Chan                  = NULL;
    delete Hit_X_start              ; Hit_X_start               = NULL;
    delete Hit_Y_start              ; Hit_Y_start               = NULL;
    delete Hit_Z_start              ; Hit_Z_start               = NULL;
    delete Hit_X_end                ; Hit_X_end                 = NULL;
    delete Hit_Y_end                ; Hit_Y_end                 = NULL;
    delete Hit_Z_end                ; Hit_Z_end                 = NULL;
    delete Hit_Time                 ; Hit_Time                  = NULL;
    delete Hit_RMS                  ; Hit_RMS                   = NULL;
    delete Hit_SADC                 ; Hit_SADC                  = NULL;
    delete Hit_Int                  ; Hit_Int                   = NULL;
    delete Hit_Peak                 ; Hit_Peak                  = NULL;
    delete Hit_True_GenType         ; Hit_True_GenType          = NULL;
    delete Hit_True_MainTrID        ; Hit_True_MainTrID         = NULL;
    delete Hit_True_EvEnergy        ; Hit_True_EvEnergy         = NULL;
    delete Hit_True_X               ; Hit_True_X                = NULL;                 
    delete Hit_True_Y               ; Hit_True_Y                = NULL;                 
    delete Hit_True_Z               ; Hit_True_Z                = NULL;                 
    delete Hit_True_Energy          ; Hit_True_Energy           = NULL;            
    delete Hit_True_nElec           ; Hit_True_nElec            = NULL;
    delete Hit_True_nIDEs           ; Hit_True_nIDEs            = NULL;
    delete Hit_True_MarleyIndex     ; Hit_True_MarleyIndex      = NULL;
                                                               
    delete Hit_AdjM5SADC            ; Hit_AdjM5SADC             = NULL;
    delete Hit_AdjM2SADC            ; Hit_AdjM2SADC             = NULL;
    delete Hit_AdjM1SADC            ; Hit_AdjM1SADC             = NULL;
    delete Hit_AdjP1SADC            ; Hit_AdjP1SADC             = NULL;
    delete Hit_AdjP2SADC            ; Hit_AdjP2SADC             = NULL;
    delete Hit_AdjP5SADC            ; Hit_AdjP5SADC             = NULL;
    delete Hit_AdjM5Chan            ; Hit_AdjM5Chan             = NULL;
    delete Hit_AdjM2Chan            ; Hit_AdjM2Chan             = NULL;
    delete Hit_AdjM1Chan            ; Hit_AdjM1Chan             = NULL;
    delete Hit_AdjP1Chan            ; Hit_AdjP1Chan             = NULL;
    delete Hit_AdjP2Chan            ; Hit_AdjP2Chan             = NULL;
    delete Hit_AdjP5Chan            ; Hit_AdjP5Chan             = NULL;
                                                               
    delete PDS_Flash_FlashID        ; PDS_Flash_FlashID         = NULL;
    delete PDS_Flash_YCenter        ; PDS_Flash_YCenter         = NULL;
    delete PDS_Flash_ZCenter        ; PDS_Flash_ZCenter         = NULL;
    delete PDS_Flash_YWidth         ; PDS_Flash_YWidth          = NULL;
    delete PDS_Flash_ZWidth         ; PDS_Flash_ZWidth          = NULL;
    delete PDS_Flash_Time           ; PDS_Flash_Time            = NULL;
    delete PDS_Flash_TimeWidth      ; PDS_Flash_TimeWidth       = NULL;
    delete PDS_Flash_TotalPE        ; PDS_Flash_TotalPE         = NULL;
    delete PDS_Flash_True_Distance  ; PDS_Flash_True_Distance   = NULL;
    delete PDS_Flash_True_GenType   ; PDS_Flash_True_GenType    = NULL;

    delete PDS_OpHit_OpChannel      ; PDS_OpHit_OpChannel       = NULL;
    delete PDS_OpHit_X              ; PDS_OpHit_X               = NULL;
    delete PDS_OpHit_Y              ; PDS_OpHit_Y               = NULL;
    delete PDS_OpHit_Z              ; PDS_OpHit_Z               = NULL;
    delete PDS_OpHit_PeakTimeAbs    ; PDS_OpHit_PeakTimeAbs     = NULL;
    delete PDS_OpHit_PeakTime       ; PDS_OpHit_PeakTime        = NULL;
    delete PDS_OpHit_Frame          ; PDS_OpHit_Frame           = NULL;
    delete PDS_OpHit_Width          ; PDS_OpHit_Width           = NULL;
    delete PDS_OpHit_Area           ; PDS_OpHit_Area            = NULL;
    delete PDS_OpHit_Amplitude      ; PDS_OpHit_Amplitude       = NULL;
    delete PDS_OpHit_PE             ; PDS_OpHit_PE              = NULL;
    delete PDS_OpHit_FastToTotal    ; PDS_OpHit_FastToTotal     = NULL;
    delete PDS_OpHit_True_GenType   ; PDS_OpHit_True_GenType    = NULL;
    delete PDS_OpHit_True_Energy    ; PDS_OpHit_True_Energy     = NULL;
    delete PDS_OpHit_True_TrackID   ; PDS_OpHit_True_TrackID    = NULL;
    delete PDS_OpHit_True_GenTypeAll; PDS_OpHit_True_GenTypeAll = NULL;
    delete PDS_OpHit_True_EnergyAll ; PDS_OpHit_True_EnergyAll  = NULL;
    delete PDS_OpHit_True_TrackIDAll; PDS_OpHit_True_TrackIDAll = NULL;
    delete PDS_OpHit_True_IndexAll  ; PDS_OpHit_True_IndexAll   = NULL;

    delete True_VertexChan          ; True_VertexChan           = NULL;
    delete True_Nu_Type             ; True_Nu_Type              = NULL;
    delete True_Nu_Lep_Type         ; True_Nu_Lep_Type          = NULL;
    delete True_Mode                ; True_Mode                 = NULL;
    delete True_CCNC                ; True_CCNC                 = NULL;
    delete True_HitNucleon          ; True_HitNucleon           = NULL;
    delete True_Target              ; True_Target               = NULL;
    delete True_MarlSample          ; True_MarlSample           = NULL;
    delete True_MarlTime            ; True_MarlTime             = NULL;
    delete True_MarlWeight          ; True_MarlWeight           = NULL;
    delete True_ENu                 ; True_ENu                  = NULL;
    delete True_ENu_Lep             ; True_ENu_Lep              = NULL;
    delete True_VertX               ; True_VertX                = NULL;
    delete True_VertY               ; True_VertY                = NULL;
    delete True_VertZ               ; True_VertZ                = NULL;
    delete True_VertexT             ; True_VertexT              = NULL;
    delete True_Px                  ; True_Px                   = NULL;
    delete True_Py                  ; True_Py                   = NULL;
    delete True_Pz                  ; True_Pz                   = NULL;
    delete True_Dirx                ; True_Dirx                 = NULL;
    delete True_Diry                ; True_Diry                 = NULL;
    delete True_Dirz                ; True_Dirz                 = NULL;
    delete True_Time                ; True_Time                 = NULL;
                                                               
    delete True_Bck_Mode            ; True_Bck_Mode             = NULL;
    delete True_Bck_VertX           ; True_Bck_VertX            = NULL;
    delete True_Bck_VertY           ; True_Bck_VertY            = NULL;
    delete True_Bck_VertZ           ; True_Bck_VertZ            = NULL;
    delete True_Bck_Time            ; True_Bck_Time             = NULL;
    delete True_Bck_Energy          ; True_Bck_Energy           = NULL;
    delete True_Bck_PDG             ; True_Bck_PDG              = NULL;
                                                               
    delete True_Bck_ModeAll         ; True_Bck_ModeAll          = NULL;
    delete True_Bck_VertXAll        ; True_Bck_VertXAll         = NULL;
    delete True_Bck_VertYAll        ; True_Bck_VertYAll         = NULL;
    delete True_Bck_VertZAll        ; True_Bck_VertZAll         = NULL;
    delete True_Bck_TimeAll         ; True_Bck_TimeAll          = NULL;
    delete True_Bck_EndVertXAll     ; True_Bck_EndVertXAll      = NULL;
    delete True_Bck_EndVertYAll     ; True_Bck_EndVertYAll      = NULL;
    delete True_Bck_EndVertZAll     ; True_Bck_EndVertZAll      = NULL;
    delete True_Bck_EndTimeAll      ; True_Bck_EndTimeAll       = NULL;
    delete True_Bck_EnergyAll       ; True_Bck_EnergyAll        = NULL;
    delete True_Bck_PDGAll          ; True_Bck_PDGAll           = NULL;
    delete True_Bck_IDAll           ; True_Bck_IDAll            = NULL;
    
    if (f_Input) f_Input->Close();

  };
  
  InputManager():
    f_Input(NULL),
    t_Input(NULL),
    Run        (0),
    SubRun     (0),
    Event      (0),
    NTotHit    (0),
    NColHit    (0),
    NIndHit    (0),
    NHitNoBT   (0),
    NFlash     (0),
    NFlashNoBT (0),
    TotGen_Marl(0),
    TotGen_APA (0),
    TotGen_CPA (0),
    TotGen_Ar39(0),
    TotGen_Neut(0),
    TotGen_Kryp(0),
    TotGen_Plon(0),
    TotGen_Rdon(0),
    TotGen_Ar42(0),   
    Hit_View                 (NULL),
    Hit_Size                 (NULL),
    Hit_TPC                  (NULL),
    Hit_Chan                 (NULL),
    Hit_X_start              (NULL),
    Hit_Y_start              (NULL),
    Hit_Z_start              (NULL),
    Hit_X_end                (NULL),
    Hit_Y_end                (NULL),
    Hit_Z_end                (NULL),
    Hit_Time                 (NULL),
    Hit_RMS                  (NULL),
    Hit_SADC                 (NULL),
    Hit_Int                  (NULL),
    Hit_Peak                 (NULL),
    Hit_True_GenType         (NULL),
    Hit_True_MainTrID        (NULL),
    Hit_True_EvEnergy        (NULL),
    Hit_True_X               (NULL),                 
    Hit_True_Y               (NULL),                 
    Hit_True_Z               (NULL),                 
    Hit_True_Energy          (NULL),            
    Hit_True_nElec           (NULL),
    Hit_True_nIDEs           (NULL),
    Hit_True_MarleyIndex     (NULL),
    
    Hit_AdjM5SADC            (NULL),
    Hit_AdjM2SADC            (NULL),
    Hit_AdjM1SADC            (NULL),
    Hit_AdjP1SADC            (NULL),
    Hit_AdjP2SADC            (NULL),
    Hit_AdjP5SADC            (NULL),
    Hit_AdjM5Chan            (NULL),
    Hit_AdjM2Chan            (NULL),
    Hit_AdjM1Chan            (NULL),
    Hit_AdjP1Chan            (NULL),
    Hit_AdjP2Chan            (NULL),
    Hit_AdjP5Chan            (NULL),

    PDS_Flash_FlashID        (NULL),
    PDS_Flash_YCenter        (NULL),
    PDS_Flash_ZCenter        (NULL),
    PDS_Flash_YWidth         (NULL),
    PDS_Flash_ZWidth         (NULL),
    PDS_Flash_Time           (NULL),
    PDS_Flash_TimeWidth      (NULL),
    PDS_Flash_TotalPE        (NULL),
    PDS_Flash_True_Distance  (NULL),
    PDS_Flash_True_GenType   (NULL),

    PDS_OpHit_OpChannel      (NULL),
    PDS_OpHit_X              (NULL),
    PDS_OpHit_Y              (NULL),
    PDS_OpHit_Z              (NULL),
    PDS_OpHit_PeakTimeAbs    (NULL),
    PDS_OpHit_PeakTime       (NULL),
    PDS_OpHit_Frame          (NULL),
    PDS_OpHit_Width          (NULL),
    PDS_OpHit_Area           (NULL),
    PDS_OpHit_Amplitude      (NULL),
    PDS_OpHit_PE             (NULL),
    PDS_OpHit_FastToTotal    (NULL),
    PDS_OpHit_True_GenType   (NULL),
    PDS_OpHit_True_Energy    (NULL),
    PDS_OpHit_True_TrackID   (NULL),
    PDS_OpHit_True_GenTypeAll(NULL),
    PDS_OpHit_True_EnergyAll (NULL),
    PDS_OpHit_True_TrackIDAll(NULL),
    PDS_OpHit_True_IndexAll  (NULL),

    True_VertexChan          (NULL),
    True_Nu_Type             (NULL),
    True_Nu_Lep_Type         (NULL),
    True_Mode                (NULL),
    True_CCNC                (NULL),
    True_HitNucleon          (NULL),
    True_Target              (NULL),
    True_MarlSample          (NULL),
    True_MarlTime            (NULL),
    True_MarlWeight          (NULL),
    True_ENu                 (NULL),
    True_ENu_Lep             (NULL),
    True_VertX               (NULL),
    True_VertY               (NULL),
    True_VertZ               (NULL),
    True_VertexT             (NULL),
    True_Px                  (NULL),
    True_Py                  (NULL),
    True_Pz                  (NULL),
    True_Dirx                (NULL),
    True_Diry                (NULL),
    True_Dirz                (NULL),
    True_Time                (NULL),

    True_Bck_Mode            (NULL),
    True_Bck_VertX           (NULL),
    True_Bck_VertY           (NULL),
    True_Bck_VertZ           (NULL),
    True_Bck_Time            (NULL),
    True_Bck_Energy          (NULL),
    True_Bck_PDG             (NULL),

    True_Bck_ModeAll         (NULL),
    True_Bck_VertXAll        (NULL),
    True_Bck_VertYAll        (NULL),
    True_Bck_VertZAll        (NULL),
    True_Bck_TimeAll         (NULL),
    True_Bck_EndVertXAll     (NULL),
    True_Bck_EndVertYAll     (NULL),
    True_Bck_EndVertZAll     (NULL),
    True_Bck_EndTimeAll      (NULL),
    True_Bck_EnergyAll       (NULL),
    True_Bck_PDGAll          (NULL),
    True_Bck_IDAll           (NULL){};
  
  int         GetEntries  ()                       const { return t_Input->GetEntries(); };
  void        GetEntry    (const int i=0)                { t_Input->GetEntry(i); };
  std::string GetInputFile()                       const { return filename; };
  std::string GetInputTree()                       const { return treename; };
  void        SetInputFile(const std::string s="")       { filename=s; };
  void        SetInputTree(const std::string s="")       { treename=s; };
    
  void LoadTree(){
    if (filename == "") {
      std::cerr << "Need to provide an input file" << std::endl;
      exit(1);
    }
    
    f_Input = new TFile(filename.c_str(), "READ");
    if (!f_Input->IsOpen()) {
      std::cerr << "The file " << filename.c_str() << " does not exist." << std::endl;
      exit(1);
    }

    if (f_Input->Get(treename.c_str())) {
      t_Input = (TTree*)f_Input->Get(treename.c_str());
    } else {
      std::cerr << "The requested tree (" << treename.c_str() << ") does not exist in file " << filename << std::endl;
      exit(1);
    }
    
    
    t_Input->SetBranchAddress("Run"       , &Run       );
    t_Input->SetBranchAddress("SubRun"    , &SubRun    );
    t_Input->SetBranchAddress("Event"     , &Event     );
    t_Input->SetBranchAddress("NTotHit"   , &NTotHit   );
    t_Input->SetBranchAddress("NColHit"   , &NColHit   );
    t_Input->SetBranchAddress("NIndHit"   , &NIndHit   );
    t_Input->SetBranchAddress("NHitNoBT"  , &NHitNoBT  );
    t_Input->SetBranchAddress("NFlash"    , &NFlash    );
    t_Input->SetBranchAddress("NFlashNoBT", &NFlashNoBT);

    t_Input->SetBranchAddress("Hit_View"                 , &Hit_View                 );
    t_Input->SetBranchAddress("Hit_Size"                 , &Hit_Size                 );
    t_Input->SetBranchAddress("Hit_TPC"                  , &Hit_TPC                  );
    t_Input->SetBranchAddress("Hit_Chan"                 , &Hit_Chan                 );
    t_Input->SetBranchAddress("Hit_X_start"              , &Hit_X_start              );
    t_Input->SetBranchAddress("Hit_Y_start"              , &Hit_Y_start              );
    t_Input->SetBranchAddress("Hit_Z_start"              , &Hit_Z_start              );
    t_Input->SetBranchAddress("Hit_X_end"                , &Hit_X_end                );
    t_Input->SetBranchAddress("Hit_Y_end"                , &Hit_Y_end                );
    t_Input->SetBranchAddress("Hit_Z_end"                , &Hit_Z_end                );
    t_Input->SetBranchAddress("Hit_Time"                 , &Hit_Time                 );
    t_Input->SetBranchAddress("Hit_RMS"                  , &Hit_RMS                  );
    t_Input->SetBranchAddress("Hit_SADC"                 , &Hit_SADC                 );
    t_Input->SetBranchAddress("Hit_Int"                  , &Hit_Int                  );
    t_Input->SetBranchAddress("Hit_Peak"                 , &Hit_Peak                 );
    t_Input->SetBranchAddress("Hit_True_GenType"         , &Hit_True_GenType         );
    t_Input->SetBranchAddress("Hit_True_MainTrID"        , &Hit_True_MainTrID        );
    t_Input->SetBranchAddress("Hit_True_EvEnergy"        , &Hit_True_EvEnergy        );
    t_Input->SetBranchAddress("Hit_True_X"               , &Hit_True_X               );
    t_Input->SetBranchAddress("Hit_True_Y"               , &Hit_True_Y               );
    t_Input->SetBranchAddress("Hit_True_Z"               , &Hit_True_Z               );
    t_Input->SetBranchAddress("Hit_True_Energy"          , &Hit_True_Energy          );
    t_Input->SetBranchAddress("Hit_True_nElec"           , &Hit_True_nElec           );
    t_Input->SetBranchAddress("Hit_True_nIDEs"           , &Hit_True_nIDEs           );
    if(t_Input->GetListOfBranches()->FindObject("Hit_True_MarleyIndex")) {
      t_Input->SetBranchAddress("Hit_True_MarleyIndex", &Hit_True_MarleyIndex);
    } else {
      std::cout << "Hit_True_MarleyIndex branch doesn't exist, assigning NULL for this vector!" << std::endl;      
      Hit_True_MarleyIndex=NULL;
    }

    // if(MyTree->GetBranch("SomeBranchName")) {
    //   MyTree->SetBranchAddress("SomeBranchName", &SomeVabiable);
    // }
    t_Input->SetBranchAddress("Hit_AdjM5SADC"            , &Hit_AdjM5SADC            );
    t_Input->SetBranchAddress("Hit_AdjM2SADC"            , &Hit_AdjM2SADC            );
    t_Input->SetBranchAddress("Hit_AdjM1SADC"            , &Hit_AdjM1SADC            );
    t_Input->SetBranchAddress("Hit_AdjP1SADC"            , &Hit_AdjP1SADC            );
    t_Input->SetBranchAddress("Hit_AdjP2SADC"            , &Hit_AdjP2SADC            );
    t_Input->SetBranchAddress("Hit_AdjP5SADC"            , &Hit_AdjP5SADC            );
    t_Input->SetBranchAddress("Hit_AdjM5Chan"            , &Hit_AdjM5Chan            );
    t_Input->SetBranchAddress("Hit_AdjM2Chan"            , &Hit_AdjM2Chan            );
    t_Input->SetBranchAddress("Hit_AdjM1Chan"            , &Hit_AdjM1Chan            );
    t_Input->SetBranchAddress("Hit_AdjP1Chan"            , &Hit_AdjP1Chan            );
    t_Input->SetBranchAddress("Hit_AdjP2Chan"            , &Hit_AdjP2Chan            );
    t_Input->SetBranchAddress("Hit_AdjP5Chan"            , &Hit_AdjP5Chan            );
                                                                                     
    t_Input->SetBranchAddress("PDS_Flash_FlashID"        , &PDS_Flash_FlashID        );
    t_Input->SetBranchAddress("PDS_Flash_YCenter"        , &PDS_Flash_YCenter        );
    t_Input->SetBranchAddress("PDS_Flash_ZCenter"        , &PDS_Flash_ZCenter        );
    t_Input->SetBranchAddress("PDS_Flash_YWidth"         , &PDS_Flash_YWidth         );
    t_Input->SetBranchAddress("PDS_Flash_ZWidth"         , &PDS_Flash_ZWidth         );
    t_Input->SetBranchAddress("PDS_Flash_Time"           , &PDS_Flash_Time           );
    t_Input->SetBranchAddress("PDS_Flash_TimeWidth"      , &PDS_Flash_TimeWidth      );
    t_Input->SetBranchAddress("PDS_Flash_TotalPE"        , &PDS_Flash_TotalPE        );
    t_Input->SetBranchAddress("PDS_Flash_True_Distance"  , &PDS_Flash_True_Distance  );
    t_Input->SetBranchAddress("PDS_Flash_True_GenType"   , &PDS_Flash_True_GenType   );
                                                                                     
    t_Input->SetBranchAddress("PDS_OpHit_OpChannel"      , &PDS_OpHit_OpChannel      );
    t_Input->SetBranchAddress("PDS_OpHit_X"              , &PDS_OpHit_X              );
    t_Input->SetBranchAddress("PDS_OpHit_Y"              , &PDS_OpHit_Y              );
    t_Input->SetBranchAddress("PDS_OpHit_Z"              , &PDS_OpHit_Z              );
    t_Input->SetBranchAddress("PDS_OpHit_PeakTimeAbs"    , &PDS_OpHit_PeakTimeAbs    );
    t_Input->SetBranchAddress("PDS_OpHit_PeakTime"       , &PDS_OpHit_PeakTime       );
    t_Input->SetBranchAddress("PDS_OpHit_Frame"          , &PDS_OpHit_Frame          );
    t_Input->SetBranchAddress("PDS_OpHit_Width"          , &PDS_OpHit_Width          );
    t_Input->SetBranchAddress("PDS_OpHit_Area"           , &PDS_OpHit_Area           );
    t_Input->SetBranchAddress("PDS_OpHit_Amplitude"      , &PDS_OpHit_Amplitude      );
    t_Input->SetBranchAddress("PDS_OpHit_PE"             , &PDS_OpHit_PE             );
    t_Input->SetBranchAddress("PDS_OpHit_FastToTotal"    , &PDS_OpHit_FastToTotal    );
    t_Input->SetBranchAddress("PDS_OpHit_True_GenType"   , &PDS_OpHit_True_GenType   );
    t_Input->SetBranchAddress("PDS_OpHit_True_Energy"    , &PDS_OpHit_True_Energy    );
    t_Input->SetBranchAddress("PDS_OpHit_True_TrackID"   , &PDS_OpHit_True_TrackID   );
    t_Input->SetBranchAddress("PDS_OpHit_True_GenTypeAll", &PDS_OpHit_True_GenTypeAll);
    t_Input->SetBranchAddress("PDS_OpHit_True_EnergyAll" , &PDS_OpHit_True_EnergyAll );
    t_Input->SetBranchAddress("PDS_OpHit_True_TrackIDAll", &PDS_OpHit_True_TrackIDAll);
    t_Input->SetBranchAddress("PDS_OpHit_True_IndexAll"  , &PDS_OpHit_True_IndexAll  );
                                                                                     
    t_Input->SetBranchAddress("True_VertexChan"          , &True_VertexChan          );
    t_Input->SetBranchAddress("True_Nu_Type"             , &True_Nu_Type             );
    t_Input->SetBranchAddress("True_Nu_Lep_Type"         , &True_Nu_Lep_Type         );
    t_Input->SetBranchAddress("True_Mode"                , &True_Mode                );
    t_Input->SetBranchAddress("True_CCNC"                , &True_CCNC                );
    t_Input->SetBranchAddress("True_HitNucleon"          , &True_HitNucleon          );
    t_Input->SetBranchAddress("True_Target"              , &True_Target              );
    t_Input->SetBranchAddress("True_MarlSample"          , &True_MarlSample          );
    t_Input->SetBranchAddress("True_MarlTime"            , &True_MarlTime            );
    t_Input->SetBranchAddress("True_MarlWeight"          , &True_MarlWeight          );
    t_Input->SetBranchAddress("True_ENu"                 , &True_ENu                 );
    t_Input->SetBranchAddress("True_ENu_Lep"             , &True_ENu_Lep             );
    t_Input->SetBranchAddress("True_VertX"               , &True_VertX               );
    t_Input->SetBranchAddress("True_VertY"               , &True_VertY               );
    t_Input->SetBranchAddress("True_VertZ"               , &True_VertZ               );
    t_Input->SetBranchAddress("True_VertexT"             , &True_VertexT             );
    t_Input->SetBranchAddress("True_Px"                  , &True_Px                  );
    t_Input->SetBranchAddress("True_Py"                  , &True_Py                  );
    t_Input->SetBranchAddress("True_Pz"                  , &True_Pz                  );
    t_Input->SetBranchAddress("True_Dirx"                , &True_Dirx                );
    t_Input->SetBranchAddress("True_Diry"                , &True_Diry                );
    t_Input->SetBranchAddress("True_Dirz"                , &True_Dirz                );
    t_Input->SetBranchAddress("True_Time"                , &True_Time                );
    
    t_Input->SetBranchAddress("True_Bck_Mode"            , &True_Bck_Mode            );
    t_Input->SetBranchAddress("True_Bck_VertX"           , &True_Bck_VertX           );
    t_Input->SetBranchAddress("True_Bck_VertY"           , &True_Bck_VertY           );
    t_Input->SetBranchAddress("True_Bck_VertZ"           , &True_Bck_VertZ           );
    t_Input->SetBranchAddress("True_Bck_Time"            , &True_Bck_Time            );
    t_Input->SetBranchAddress("True_Bck_Energy"          , &True_Bck_Energy          );
    t_Input->SetBranchAddress("True_Bck_PDG"             , &True_Bck_PDG             );

    t_Input->SetBranchAddress("True_Bck_ModeAll"         , &True_Bck_ModeAll         );
    t_Input->SetBranchAddress("True_Bck_VertXAll"        , &True_Bck_VertXAll        );
    t_Input->SetBranchAddress("True_Bck_VertYAll"        , &True_Bck_VertYAll        );
    t_Input->SetBranchAddress("True_Bck_VertZAll"        , &True_Bck_VertZAll        );
    t_Input->SetBranchAddress("True_Bck_TimeAll"         , &True_Bck_TimeAll         );
    t_Input->SetBranchAddress("True_Bck_EndVertXAll"     , &True_Bck_EndVertXAll     );
    t_Input->SetBranchAddress("True_Bck_EndVertYAll"     , &True_Bck_EndVertYAll     );
    t_Input->SetBranchAddress("True_Bck_EndVertZAll"     , &True_Bck_EndVertZAll     );
    t_Input->SetBranchAddress("True_Bck_EndTimeAll"      , &True_Bck_EndTimeAll      );
    t_Input->SetBranchAddress("True_Bck_EnergyAll"       , &True_Bck_EnergyAll       );
    t_Input->SetBranchAddress("True_Bck_PDGAll"          , &True_Bck_PDGAll          );
    t_Input->SetBranchAddress("True_Bck_IDAll"           , &True_Bck_IDAll           );
    
    t_Input->SetBranchAddress("TotGen_Marl", &TotGen_Marl);
    t_Input->SetBranchAddress("TotGen_APA" , &TotGen_APA );
    t_Input->SetBranchAddress("TotGen_CPA" , &TotGen_CPA );
    t_Input->SetBranchAddress("TotGen_Ar39", &TotGen_Ar39);
    t_Input->SetBranchAddress("TotGen_Neut", &TotGen_Neut);
    t_Input->SetBranchAddress("TotGen_Kryp", &TotGen_Kryp);
    t_Input->SetBranchAddress("TotGen_Plon", &TotGen_Plon);
    t_Input->SetBranchAddress("TotGen_Rdon", &TotGen_Rdon);
    t_Input->SetBranchAddress("TotGen_Ar42", &TotGen_Ar42);
    
  };
  
  
};
  
#endif

