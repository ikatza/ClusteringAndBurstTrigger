#ifndef CLUSTERENERGYESTIMATOR_HH
#define CLUSTERENERGYESTIMATOR_HH
#include "WireCluster.hh"

#include "TXMLEngine.h"
#include "TMath.h"

#include <stdio.h>
#include <stdlib.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

class ClusterEnergyEstimator {

public:
  ClusterEnergyEstimator(const std::string filename);
  ~ClusterEnergyEstimator();
  double EstimateEnergy(const WireCluster&) const;
  void   EstimateEnergy(const std::vector<WireCluster*>&) const;
  void   EstimateEnergy(const std::vector<WireCluster*>*) const;
  std::vector<double> GetElapsedTime() const { return {0}; };
  enum Method{
    kUnknown   =-1,
    kLinReg    = 0,
    kNeuralNet = 1
  };
  
private:
  void ParseXMLFile();
  void ParseXMLFileLinearRegression(std::map<std::string,XMLNodePointer_t>);
  void ParseXMLFileNeuralNetwork(std::map<std::string,XMLNodePointer_t>);

  void ActivationFunctionRelu(boost::numeric::ublas::vector<double>& v) const {
    for (size_t i=0; i<v.size(); ++i) {
      v[i] = TMath::Max(0.,v[i]);
    }
  };
  void ActivationFunctionSigmoid(boost::numeric::ublas::vector<double>& v) const {
    for (size_t i=0; i<v.size(); ++i) {
      v[i] = 1. / (1. + TMath::Exp(-v[i]));
    }
  };
  
  void LoadChildAndAssertExistance(const XMLNodePointer_t a,
                                   std::map<std::string,XMLNodePointer_t>& b,
                                   const std::vector<std::string>& c,
                                   const int i=1) const;
  boost::numeric::ublas::vector<double> NeuralNetwork   (const boost::numeric::ublas::vector<double>&) const;
  boost::numeric::ublas::vector<double> LinearRegression(const boost::numeric::ublas::vector<double>&) const;
  void PrintTestCluster() const;

  void NormaliseFeature(boost::numeric::ublas::vector<double>& ,
                        const boost::numeric::ublas::vector<double>,
                        const boost::numeric::ublas::vector<double>) const;
  void RecoverFeature(boost::numeric::ublas::vector<double>&,
                      const boost::numeric::ublas::vector<double>,
                      const boost::numeric::ublas::vector<double>) const;
  TXMLEngine* fXML;
  std::string fFileName;
  int fNNeuron;
  int fNParam;
  int fNOutput;
  boost::numeric::ublas::matrix<double> fWeight_h1        ;
  boost::numeric::ublas::matrix<double> fWeight_out       ;
  boost::numeric::ublas::matrix<double> fBias_b1          ;
  boost::numeric::ublas::matrix<double> fBias_out         ;
  boost::numeric::ublas::vector<double> fInputData_Mean   ;
  boost::numeric::ublas::vector<double> fInputData_StdDev ;
  boost::numeric::ublas::vector<double> fOutputData_Mean  ;
  boost::numeric::ublas::vector<double> fOutputData_StdDev;
  Method fMethod;
  std::vector<double> fTimeElapsed;
  
};

#endif
