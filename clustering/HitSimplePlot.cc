#include "HitSimplePlot.hh"

HitSimplePlot::HitSimplePlot(){
  im.SetInputFile("/dune/app/users/plasorak/workspace/SNAna.root");
  im.SetInputTree("snanagaushit/SNSimTree");
  im.LoadTree();
  c = new TCanvas;
  nbin_1D_RMS    [1] = 100 ; nbin_1D_RMS    [2] = 40  ; nbin_1D_RMS    [3] = 100 ; nbin_1D_RMS    [4] = 20 ;
  xmin_1D_RMS    [1] = 0   ; xmin_1D_RMS    [2] = 0   ; xmin_1D_RMS    [3] = 0   ; xmin_1D_RMS    [4] = 0   ;
  xmax_1D_RMS    [1] = 4000; xmax_1D_RMS    [2] = 200 ; xmax_1D_RMS    [3] = 1000; xmax_1D_RMS    [4] = 20;
  nbin_1D_Width  [1] = 100 ; nbin_1D_Width  [2] = 40  ; nbin_1D_Width  [3] = 100 ; nbin_1D_Width  [4] = 20 ;
  xmin_1D_Width  [1] = 0   ; xmin_1D_Width  [2] = 0   ; xmin_1D_Width  [3] = 0   ; xmin_1D_Width  [4] = 0   ;
  xmax_1D_Width  [1] = 4000; xmax_1D_Width  [2] = 1000; xmax_1D_Width  [3] = 1000; xmax_1D_Width  [4] = 20;
  nbin_2D_x_RMS  [1] = 100 ; nbin_2D_x_RMS  [2] = 40  ; nbin_2D_x_RMS  [3] = 100 ; nbin_2D_x_RMS  [4] = 20 ;
  xmin_2D_x_RMS  [1] = 0   ; xmin_2D_x_RMS  [2] = 0   ; xmin_2D_x_RMS  [3] = 0   ; xmin_2D_x_RMS  [4] = 0   ;
  xmax_2D_x_RMS  [1] = 4000; xmax_2D_x_RMS  [2] = 200 ; xmax_2D_x_RMS  [3] = 1000; xmax_2D_x_RMS  [4] = 20;
  nbin_2D_x_Width[1] = 100 ; nbin_2D_x_Width[2] = 40  ; nbin_2D_x_Width[3] = 100 ; nbin_2D_x_Width[4] = 20 ;
  xmin_2D_x_Width[1] = 0   ; xmin_2D_x_Width[2] = 0   ; xmin_2D_x_Width[3] = 0   ; xmin_2D_x_Width[4] = 0   ;
  xmax_2D_x_Width[1] = 4000; xmax_2D_x_Width[2] = 1000; xmax_2D_x_Width[3] = 1000; xmax_2D_x_Width[4] = 20;
  nbin_2D_y_RMS  [1] = 100 ; nbin_2D_y_RMS  [2] = 40  ; nbin_2D_y_RMS  [3] = 100 ; nbin_2D_y_RMS  [4] = 20 ;
  xmin_2D_y_RMS  [1] = 0   ; xmin_2D_y_RMS  [2] = 0   ; xmin_2D_y_RMS  [3] = 0   ; xmin_2D_y_RMS  [4] = 0   ;
  xmax_2D_y_RMS  [1] = 4000; xmax_2D_y_RMS  [2] = 200 ; xmax_2D_y_RMS  [3] = 1000; xmax_2D_y_RMS  [4] = 20;
  nbin_2D_y_Width[1] = 100 ; nbin_2D_y_Width[2] = 40  ; nbin_2D_y_Width[3] = 100 ; nbin_2D_y_Width[4] = 20 ;
  xmin_2D_y_Width[1] = 0   ; xmin_2D_y_Width[2] = 0   ; xmin_2D_y_Width[3] = 0   ; xmin_2D_y_Width[4] = 0   ;
  xmax_2D_y_Width[1] = 2000; xmax_2D_y_Width[2] = 1000; xmax_2D_y_Width[3] = 1000; xmax_2D_y_Width[4] = 20;
};


HitSimplePlot::~HitSimplePlot(){
  delete c;
  c = NULL;
};


void HitSimplePlot::GetDimension(const int inDim, int& dim1, int& dim2) const {
  if (inDim<0) {
    std::cerr << "Invalid dimension!" << std::endl;
    exit(1);
  }
  
  if (inDim<10) {
    dim1 = inDim;
    dim2 = -1;
  }

  if (inDim>=10 && inDim<20) {
    dim1 = 1;
    dim2 = inDim-10;
  }

  if (inDim>=20 && inDim<30) {
    dim1 = 2;
    dim2 = inDim-20;
  }

  if (inDim>=30 && inDim<40) {
    dim1 = 3;
    dim2 = inDim-30;
  }

  if (inDim>=40 && inDim<50) {
    dim1 = 4;
    dim2 = inDim-40;
  }
};


void HitSimplePlot::PlotRMS(const int dim) {
  int dim1, dim2;
  GetDimension(dim, dim1, dim2);
  if (dim2 == -1) { PlotRMS1D(dim1); }
  else            { PlotRMS2D(dim1,dim2); }
};


void HitSimplePlot::PlotWidth(const int dim) {
  int dim1, dim2;
  GetDimension(dim, dim1, dim2);
  if (dim2 == -1) { PlotWidth1D(dim1); }
  else            { PlotWidth2D(dim1,dim2); }
};


double HitSimplePlot::GetMean(std::vector<double> const vec) const {
  double mean=0;
  for (auto const& it : vec) {
    mean+=it;
  }
  mean/=(double)vec.size();
  return mean;
};


double HitSimplePlot::GetStdDev(std::vector<double> const vec) const {
  double stddev=0;
  double mean=GetMean(vec);
  for (auto const& it : vec)
    stddev+=(mean-it)*(mean-it);
  stddev/=(double)vec.size();
  stddev = sqrt(stddev);
  return stddev;
};


void HitSimplePlot::GetMinMax(std::vector<double> const vec, double& min, double& max) const {
  min=*std::min_element(vec.begin(), vec.end());
  max=*std::max_element(vec.begin(), vec.end());
  return;
};


std::map<GenType, std::map<int,std::vector<double> > > HitSimplePlot::GetMap(std::vector<double> const* pos,
                                                                             std::vector<int> const* gen,
                                                                             std::vector<int> const* ids) const {
  std::map<GenType, std::map<int, std::vector<double> > > map_vec_pos;
  for (size_t it=0; it<gen->size(); ++it) {
    GenType gent = ConvertIntToGenType((*gen)[it]);
    if (gent == kSNnu) { map_vec_pos[gent][0].push_back((*pos)[it]); }
    else               { map_vec_pos[gent][(*ids)[it]].push_back((*pos)[it]); }
  }
 
  return map_vec_pos;
};


std::map<GenType, std::map<int,std::pair<std::vector<double>,std::vector<double> > > > HitSimplePlot::Get2DMap(
  std::vector<double> const* pos1,
  std::vector<double> const* pos2,
  std::vector<int> const* gen,
  std::vector<int> const* ids) const {
  std::map<GenType, std::map<int,std::pair<std::vector<double>,std::vector<double> > > >map_vec_pos;
  for (size_t it=0; it<gen->size(); ++it) {
    GenType gent = ConvertIntToGenType((*gen)[it]);
    if (gent == kSNnu) {
      map_vec_pos[gent][0].first .push_back((*pos1)[it]);
      map_vec_pos[gent][0].second.push_back((*pos2)[it]);
    } else {
      map_vec_pos[gent][(*ids)[it]].first .push_back((*pos1)[it]);
      map_vec_pos[gent][(*ids)[it]].second.push_back((*pos2)[it]);
    }
  }
 
  return map_vec_pos;
};


void HitSimplePlot::FillHistosRMS(std::vector<double>* vector, std::map<GenType, TH1D*>& RMS) {
  for (int iEntry=0; iEntry<im.GetEntries(); ++iEntry) {
    PrintProgress(iEntry,im.GetEntries());
    im.GetEntry(iEntry);
    std::map<GenType,std::map<int, std::vector<double> > > map_pos = GetMap(vector,
                                                                            im.PDS_OpHit_True_GenType,
                                                                            im.PDS_OpHit_True_TrackID);
    for (auto const& it: map_pos) {
      for (auto const& it2: it.second) {
        double stddev = GetStdDev(it2.second);
        if(it2.second.size()>1) RMS[it.first]->Fill(stddev);
      }
    }
  }
};


void HitSimplePlot::Fill2DHistosRMS(std::vector<double>* vecx,
                                    std::vector<double>* vecy,
                                    std::map<GenType, TH2D*>& RMS) {
  
  for (int iEntry=0; iEntry<im.GetEntries(); ++iEntry) {
    PrintProgress(iEntry,im.GetEntries());
    im.GetEntry(iEntry);
    std::map<GenType,std::map<int, std::pair<std::vector<double>,std::vector<double> > > > map_pos =
      Get2DMap(vecx, vecy,
               im.PDS_OpHit_True_GenType,
               im.PDS_OpHit_True_TrackID);
    
    for (auto const& it: map_pos) {
      for (auto const& it2: it.second) {
        double stddev1 = GetStdDev(it2.second.first);
        double stddev2 = GetStdDev(it2.second.second);
        if(it2.second.second.size()>1) RMS[it.first]->Fill(stddev1, stddev2);
      }
    }
  }
};


void HitSimplePlot::FillHistosWidth(std::vector<double>* vector, std::map<GenType, TH1D*>& Width) {
  for (int iEntry=0; iEntry<im.GetEntries(); ++iEntry) {
    PrintProgress(iEntry,im.GetEntries());
    im.GetEntry(iEntry);
    std::map<GenType,std::map<int, std::vector<double> > > map_pos = GetMap(vector,
                                                                            im.PDS_OpHit_True_GenType,
                                                                            im.PDS_OpHit_True_TrackID);
    for (auto const& it: map_pos) {
      for (auto const& it2: it.second) {
        double min;
        double max;
        GetMinMax(it2.second, min, max);
        if(it2.second.size()>1) Width[it.first]->Fill(abs(max-min));
      }
    }
  }
};


void HitSimplePlot::Fill2DHistosWidth(std::vector<double>* vecx,
                                      std::vector<double>* vecy,
                                      std::map<GenType, TH2D*>& Width) {
  for (int iEntry=0; iEntry<im.GetEntries(); ++iEntry) {
    PrintProgress(iEntry,im.GetEntries());
    im.GetEntry(iEntry);
    std::map<GenType,std::map<int, std::pair<std::vector<double>,std::vector<double> > > > map_pos =
      Get2DMap(vecx, vecy,
               im.PDS_OpHit_True_GenType,
               im.PDS_OpHit_True_TrackID);
    for (auto const& it: map_pos) {
      for (auto const& it2: it.second) {
        double min1;
        double max1;
        double min2;
        double max2;
        GetMinMax(it2.second.first , min1, max1);
        GetMinMax(it2.second.second, min2, max2);
        if(it2.second.second.size()>1) Width[it.first]->Fill(abs(max1-min1),abs(max2-min2));
      }
    }
  }
};


void HitSimplePlot::PlotRMS1D(const int dim) {
  std::string Title="";

  if      (dim==1) { Title = ";X RMS [cm];Event"; }
  else if (dim==2) { Title = ";Y RMS [cm];Event"; }
  else if (dim==3) { Title = ";Z RMS [cm];Event"; }
  else if (dim==4) { Title = ";Time RMS [#mus];Event"; }
  else {
    std::cerr << "Incorrect dim" << std::endl;
    exit(1);
  }
  
  std::map<GenType, TH1D*> RMS = GetHistos(Form("RMS_d%i",dim),
                                           Title,
                                           nbin_1D_RMS[dim],xmin_1D_RMS[dim],xmax_1D_RMS[dim]);

  if      (dim==1) { FillHistosRMS(im.PDS_OpHit_X,        RMS); }
  else if (dim==2) { FillHistosRMS(im.PDS_OpHit_Y,        RMS); }
  else if (dim==3) { FillHistosRMS(im.PDS_OpHit_Z,        RMS); }
  else if (dim==4) { FillHistosRMS(im.PDS_OpHit_PeakTime, RMS); }

  c->cd();
  std::string filename = Form("RMS_1D_%i.pdf",dim);
  c->Print((filename+"[").c_str());
  for (auto const& it: RMS) {
    AddOverflow(it.second);
    it.second->Draw();
    c->Print(filename.c_str());
  }
  c->Print((filename+"]").c_str());
           
};


void HitSimplePlot::PlotWidth1D(const int dim) {
   std::string Title="";

  if      (dim==1) { Title = ";X Width [cm];Event"; }
  else if (dim==2) { Title = ";Y Width [cm];Event"; }
  else if (dim==3) { Title = ";Z Width [cm];Event"; }
  else if (dim==4) { Title = ";Time Width [#mus];Event"; }
  else {
    std::cerr << "Incorrect dim" << std::endl;
    exit(1);
  }
  
  std::map<GenType, TH1D*> Width = GetHistos(Form("Width_d%i",dim),
                                             Title,
                                             nbin_1D_Width[dim],xmin_1D_Width[dim],xmax_1D_Width[dim]);

  if      (dim==1) { FillHistosWidth(im.PDS_OpHit_X,        Width); }
  else if (dim==2) { FillHistosWidth(im.PDS_OpHit_Y,        Width); }
  else if (dim==3) { FillHistosWidth(im.PDS_OpHit_Z,        Width); }
  else if (dim==4) { FillHistosWidth(im.PDS_OpHit_PeakTime, Width); }

  c->cd();
  std::string filename = Form("Width_1D_%i.pdf",dim);
  c->Print((filename+"[").c_str());
  for (auto const& it: Width) {
    AddOverflow(it.second);
    it.second->Draw();
    c->Print(filename.c_str());
  }
  c->Print((filename+"]").c_str());
};


void HitSimplePlot::PlotRMS2D(const int dim1, const int dim2) {
  std::string Title="";

  std::vector<double>* vecx;
  if      (dim1==1) { Title = ";X RMS [cm]";      vecx = im.PDS_OpHit_X;       }
  else if (dim1==2) { Title = ";Y RMS [cm]";      vecx = im.PDS_OpHit_Y;       }
  else if (dim1==3) { Title = ";Z RMS [cm]";      vecx = im.PDS_OpHit_Z;       }
  else if (dim1==4) { Title = ";Time RMS [#mus]"; vecx = im.PDS_OpHit_PeakTime;}
  else {
    std::cerr << "Incorrect dim" << std::endl;
    exit(1);
  }
  
  std::vector<double>* vecy;
  if      (dim2==1) { Title += ";X RMS [cm]";      vecy = im.PDS_OpHit_X;       }
  else if (dim2==2) { Title += ";Y RMS [cm]";      vecy = im.PDS_OpHit_Y;       }
  else if (dim2==3) { Title += ";Z RMS [cm]";      vecy = im.PDS_OpHit_Z;       }
  else if (dim2==4) { Title += ";Time RMS [#mus]"; vecy = im.PDS_OpHit_PeakTime;}
  else {
    std::cerr << "Incorrect dim" << std::endl;
    exit(1);
  }
  
  std::map<GenType, TH2D*> RMS = Get2DHistos(Form("RMS_d%i_%i",dim1,dim2),
                                             Title,
                                             nbin_2D_x_RMS[dim1],xmin_2D_x_RMS[dim1],xmax_2D_x_RMS[dim1],
                                             nbin_2D_y_RMS[dim2],xmin_2D_y_RMS[dim2],xmax_2D_y_RMS[dim2]);

  Fill2DHistosRMS(vecx,vecy,RMS);

  c->cd();
  std::string filename = Form("RMS_2D_%i_%i.pdf",dim1,dim2);
  c->Print((filename+"[").c_str());
  for (auto const& it: RMS) {
    it.second->Draw("COLZ");
    c->Print(filename.c_str());
  }
  c->Print((filename+"]").c_str());
  
};

void HitSimplePlot::PlotWidth2D(const int dim1, const int dim2) {
  std::string Title="";

  std::vector<double>* vecx;
  if      (dim1==1) { Title = ";X Width [cm]";      vecx = im.PDS_OpHit_X;       }
  else if (dim1==2) { Title = ";Y Width [cm]";      vecx = im.PDS_OpHit_Y;       }
  else if (dim1==3) { Title = ";Z Width [cm]";      vecx = im.PDS_OpHit_Z;       }
  else if (dim1==4) { Title = ";Time Width [#mus]"; vecx = im.PDS_OpHit_PeakTime;}
  else {
    std::cerr << "Incorrect dim" << std::endl;
    exit(1);
  }
  
  std::vector<double>* vecy;
  if      (dim2==1) { Title += ";X Width [cm]";      vecy = im.PDS_OpHit_X;       }
  else if (dim2==2) { Title += ";Y Width [cm]";      vecy = im.PDS_OpHit_Y;       }
  else if (dim2==3) { Title += ";Z Width [cm]";      vecy = im.PDS_OpHit_Z;       }
  else if (dim2==4) { Title += ";Time Width [#mus]"; vecy = im.PDS_OpHit_PeakTime;}
  else {
    std::cerr << "Incorrect dim" << std::endl;
    exit(1);
  }
  
  std::map<GenType, TH2D*> Width = Get2DHistos(Form("Width_d%i_%i",dim1,dim2),
                                               Title,
                                               nbin_2D_x_Width[dim1],xmin_2D_x_Width[dim1],xmax_2D_x_Width[dim1],
                                               nbin_2D_y_Width[dim2],xmin_2D_y_Width[dim2],xmax_2D_y_Width[dim2]);

  Fill2DHistosWidth(vecx,vecy,Width);

  c->cd();
  std::string filename = Form("Width_2D_%i_%i.pdf",dim1,dim2);
  c->Print((filename+"[").c_str());
  for (auto const& it: Width) {
    it.second->Draw("COLZ");
    c->Print(filename.c_str());
  }
  c->Print((filename+"]").c_str());
  
};
