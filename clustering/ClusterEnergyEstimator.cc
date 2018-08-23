#include "ClusterEnergyEstimator.hh"

using namespace boost::numeric::ublas;

ClusterEnergyEstimator::ClusterEnergyEstimator(std::string f):
  fXML(NULL),
  fFileName(f),
  fNNeuron(0),
  fNParam(0),
  fNOutput(0),
  fWeight_h1 (),
  fWeight_out(),
  fBias_b1   (),
  fBias_out  (),
  fInputData_Mean   (),
  fInputData_StdDev (),
  fOutputData_Mean  (),
  fOutputData_StdDev(),
  fMethod(kUnknown),
  fTimeElapsed()
{
  ParseXMLFile();
};

ClusterEnergyEstimator::~ClusterEnergyEstimator() {
  if (fXML) delete fXML;
  fXML=NULL;
  fWeight_h1 .clear();
  fWeight_out.clear();
  fBias_b1   .clear();
  fBias_out  .clear();
  fInputData_Mean   .clear();
  fInputData_StdDev .clear();
  fOutputData_Mean  .clear();
  fOutputData_StdDev.clear();
  fTimeElapsed.clear();
};

double ClusterEnergyEstimator::EstimateEnergy(const WireCluster& cluster) const {

  vector<double> Feature(fNParam);
  if (fNParam == 1) {
    Feature[0] = cluster.GetHitSumADC();
  } else {
    Feature[0] = static_cast<double>(cluster.GetChannelWidth());
    Feature[1] = static_cast<double>(cluster.GetNHit());
    Feature[2] = cluster.GetHitSumADC();
    Feature[3] = cluster.GetTimeWidth();
  }
  NormaliseFeature(Feature, fInputData_Mean, fInputData_StdDev);
  vector<double> ereco(1);
  switch (fMethod) {
  case kNeuralNet:
    ereco = NeuralNetwork(Feature);
    break;
  case kLinReg:
    ereco = LinearRegression(Feature);
    break;
  case kUnknown:
    std::cout << "Unknown energy estimation method!!" << std::endl;
    exit(1);
  }
  RecoverFeature(ereco, fOutputData_Mean, fOutputData_StdDev);
  
  return ereco(0);
};

vector<double> ClusterEnergyEstimator::NeuralNetwork(const vector<double>& Feature) const {

  vector<double> layer_1(fNNeuron);
  layer_1 = prod(Feature,fWeight_h1);
  layer_1 += column(fBias_b1,0);
  ActivationFunctionRelu(layer_1);
  vector<double> output_layer = prod(layer_1,fWeight_out);
  output_layer += row(fBias_out,0);
  return output_layer;
};

vector<double> ClusterEnergyEstimator::LinearRegression(const vector<double>& Feature) const {

  vector<double> output(1);
  output = prod(Feature,fWeight_out);
  output += row(fBias_out,0);
  return output;
};


void ClusterEnergyEstimator::NormaliseFeature(vector<double>& Feature,
                                              const vector<double> Mean,
                                              const vector<double> StdDev) const {
  Feature = element_div(Feature-Mean, StdDev);
};

void ClusterEnergyEstimator::RecoverFeature(vector<double>& Feature,
                                            const vector<double> Mean,
                                            const vector<double> StdDev) const {
  Feature = element_prod(Feature, StdDev) + Mean;
};


void ClusterEnergyEstimator::EstimateEnergy(const std::vector<WireCluster*>& clusters) const {
  for (auto& it : clusters)
    it->SetEReco(EstimateEnergy(*it));
  return;
};

void ClusterEnergyEstimator::EstimateEnergy(const std::vector<WireCluster*>* clusters) const{
  for (auto& it : (*clusters))
    it->SetEReco(EstimateEnergy(*it));
  return;
};

void ClusterEnergyEstimator::ParseXMLFile() {

  std::cout << "Parsing the XML file" << std::endl;
  // First create engine
  fXML = new TXMLEngine();

  // Now try to parse xml file
  // Only file with restricted xml syntax are supported
  XMLDocPointer_t xmldoc = fXML->ParseFile(fFileName.c_str());
  if (!xmldoc) return;
  
  // take access to main node
  XMLNodePointer_t mainnode = fXML->DocGetRootElement(xmldoc);

  std::map<std::string,XMLNodePointer_t> map_main;
  std::vector<std::string> requiredchild = {"Param","InputData","OutputData","Bias","Weight"};
  LoadChildAndAssertExistance(mainnode, map_main, requiredchild);
  std::map<std::string,XMLNodePointer_t> map_param;
  requiredchild = {"nNeuron","nParam","Activation","nOutput"};
  LoadChildAndAssertExistance(map_main["Param"], map_param, requiredchild);
  std::string content = fXML->GetNodeContent(map_param.at("nNeuron"));
  fNNeuron = stoi(content);
  content = fXML->GetNodeContent(map_param.at("nOutput"));
  fNOutput = stoi(content);
  if (fNOutput > 1) {
    std::cout << "Doesnt support more than 1 output!" << std::endl;
    exit(1);
  }
  content = fXML->GetNodeContent(map_param.at("nParam"));
  fNParam = stoi(content);
  
  std::map<std::string,XMLNodePointer_t> map_inputdata;
  std::map<std::string,XMLNodePointer_t> map_outputdata;

  if (fNParam > 1) {
    requiredchild = {"Mean","StdDev"};
    LoadChildAndAssertExistance(map_main["InputData"] , map_inputdata,  requiredchild, fNParam);
  } else {
    requiredchild = {"Mean_0","StdDev_0"};
    LoadChildAndAssertExistance(map_main["InputData"] , map_inputdata,  requiredchild, fNParam);
  }
  requiredchild = {"Mean","StdDev"};
  LoadChildAndAssertExistance(map_main["OutputData"], map_outputdata, requiredchild, fNOutput);
  fInputData_Mean    = vector<double>(fNParam);
  fInputData_StdDev  = vector<double>(fNParam);
  fOutputData_Mean   = vector<double>(fNOutput);
  fOutputData_StdDev = vector<double>(fNOutput);
  
  for (int param=0; param<fNParam; ++param) {
    std::string meanstr= Form("Mean_%i",param);
    if (map_inputdata.find(meanstr) != map_inputdata.end()) {
      content = fXML->GetNodeContent(map_inputdata.at(meanstr));
      double d = stod(content);
      fInputData_Mean(param) = d;
    }

    std::string stddevstr= Form("StdDev_%i",param);
    if (map_inputdata.find(stddevstr) != map_inputdata.end()) {
      content = fXML->GetNodeContent(map_inputdata.at(stddevstr));
      double d = stod(content);
      fInputData_StdDev(param) = d;
    }
  }

  std::string meanstr= "Mean";
  if (map_outputdata.find(meanstr) != map_outputdata.end()) {
    content = fXML->GetNodeContent(map_outputdata.at(meanstr));
    double d = stod(content);
    fOutputData_Mean(0) = d;
  }

  std::string stddevstr= "StdDev";
  if (map_outputdata.find(stddevstr) != map_outputdata.end()) {
    content = fXML->GetNodeContent(map_outputdata.at(stddevstr));
    double d = stod(content);
    fOutputData_StdDev(0) = d;
  }
  
  
  if (fNNeuron == fNParam) {
    std::cout << "Using linear regression with " << fNParam << " variables!" << std::endl;
    fMethod = kLinReg;
    ParseXMLFileLinearRegression(map_main);
  } else {
    std::cout << "Using neural network with " << fNParam << " variables and " << fNNeuron << " neurons." << std::endl;
    fMethod = kNeuralNet;
    ParseXMLFileNeuralNetwork(map_main);
  }
  
  fXML->FreeDoc(xmldoc);
  PrintTestCluster();

};



void ClusterEnergyEstimator::PrintTestCluster() const {

  switch (fMethod) {
  case kUnknown:
    std::cout << "Unknown energy estimation method!!" << std::endl;
    exit(1);
  case kLinReg:
    std::cout << "EReco on test cluster not implemented yet!!" << std::endl;
    break;
  case kNeuralNet:
    vector<double> Feature(fNParam);
    Feature[0] = 8;
    Feature[1] = 7;
    Feature[2] = 1000;
    Feature[3] = 1;
  
    NormaliseFeature(Feature, fInputData_Mean, fInputData_StdDev);
    vector<double> ereco(1);
    ereco = NeuralNetwork(Feature);
    RecoverFeature(ereco, fOutputData_Mean, fOutputData_StdDev);
    std::cout << "EReco on test cluster: "  << ereco(0) << std::endl;
    break;
  }
  
};


void ClusterEnergyEstimator::ParseXMLFileLinearRegression(std::map<std::string,XMLNodePointer_t> map_main) {
  std::string content;

  std::map<std::string,XMLNodePointer_t> map_bias;
  std::vector<std::string> requiredchild = {"out"};
  LoadChildAndAssertExistance(map_main["Bias"], map_bias, requiredchild);

  std::map<std::string,XMLNodePointer_t> map_bias_out;
  requiredchild = {"Layer_1_0"};
  LoadChildAndAssertExistance(map_bias["out"], map_bias_out, requiredchild);

  std::map<std::string,XMLNodePointer_t> map_weight;
  requiredchild = {"out"};
  LoadChildAndAssertExistance(map_main["Weight"], map_weight, requiredchild);

  std::map<std::string,XMLNodePointer_t> map_weight_out_param;
  if (fNParam==1) {
    requiredchild = {"Param_0"};
    LoadChildAndAssertExistance(map_weight["out"], map_weight_out_param, requiredchild);
  } else {
    requiredchild = {"Param"};
    LoadChildAndAssertExistance(map_weight["out"], map_weight_out_param, requiredchild, fNParam);
  }
  std::map<std::string,std::map<std::string,XMLNodePointer_t>> map_weight_out_param_layer;
  for(int i=0; i<fNParam; ++i) {
    requiredchild = {"Layer_1_0"};
    std::map<std::string,XMLNodePointer_t> dummy;
    LoadChildAndAssertExistance(map_weight_out_param.at(Form("Param_%i",i)), dummy, requiredchild);
    map_weight_out_param_layer[Form("Param_%i",i)] = dummy;
  }

  fWeight_out = matrix<double>(fNParam,fNOutput);
  fBias_out   = matrix<double>(1,fNOutput);

  for (int i=0; i<fNParam; ++i) {
    for (int j=0; j<fNOutput; ++j) {
      std::string paramstr= Form("Param_%i",i);
    
      if (map_weight_out_param_layer.find(paramstr) != map_weight_out_param_layer.end()) {
        if (map_weight_out_param_layer.at(paramstr).find("Layer_1_0") != map_weight_out_param_layer.at(paramstr).end()) {
          content = fXML->GetNodeContent(map_weight_out_param_layer.at(paramstr).at("Layer_1_0"));
          double d = stod(content);
          fWeight_out(i,j) = d;
        }
      }
    }
  }
  
  content = fXML->GetNodeContent(map_bias_out["Layer_1_0"]);
  double d = stod(content);
  fBias_out(0,0) = d;
};


void ClusterEnergyEstimator::ParseXMLFileNeuralNetwork(std::map<std::string,XMLNodePointer_t> map_main) {

  std::map<std::string,XMLNodePointer_t> map_bias;
  std::vector<std::string> requiredchild = {"b1","out"};
  std::string content;
  LoadChildAndAssertExistance(map_main["Bias"], map_bias, requiredchild);

  std::map<std::string,XMLNodePointer_t> map_bias_b1;
  requiredchild = {"Layer_1"};
  LoadChildAndAssertExistance(map_bias["b1"], map_bias_b1, requiredchild, fNNeuron);

  std::map<std::string,XMLNodePointer_t> map_bias_out;
  requiredchild = {"Layer_1_0"};
  LoadChildAndAssertExistance(map_bias["out"], map_bias_out, requiredchild);

  std::map<std::string,XMLNodePointer_t> map_weight;
  requiredchild = {"h1","out"};
  LoadChildAndAssertExistance(map_main["Weight"], map_weight, requiredchild);

  std::map<std::string,XMLNodePointer_t> map_weight_h1_param;
  requiredchild = {"Param"};
  LoadChildAndAssertExistance(map_weight["h1"], map_weight_h1_param, requiredchild, fNParam);
  
  std::map<std::string,std::map<std::string,XMLNodePointer_t>> map_weight_h1_param_layer;
  requiredchild = {"Layer_1"};
  for(int i=0; i<fNParam; ++i) {
    std::map<std::string,XMLNodePointer_t> dummy;
    LoadChildAndAssertExistance(map_weight_h1_param.at(Form("Param_%i",i)),
                                dummy, requiredchild, fNNeuron);
    map_weight_h1_param_layer[Form("Param_%i",i)] = dummy;
  }
  
  std::map<std::string,XMLNodePointer_t> map_weight_out;
  requiredchild = {"Layer_1"};
  LoadChildAndAssertExistance(map_weight["out"], map_weight_out, requiredchild, fNNeuron);

  fWeight_h1  = matrix<double>(fNParam,fNNeuron);
  fWeight_out = matrix<double>(fNNeuron,1);
  fBias_b1    = matrix<double>(fNNeuron,1);
  fBias_out   = matrix<double>(fNOutput,1);
  
  for (int param=0; param<fNParam; ++param) {
    std::string paramstr= Form("Param_%i",param);
    
    if (map_weight_h1_param_layer.find(paramstr) != map_weight_h1_param_layer.end()) {
      
      for (int neur=0; neur<fNNeuron; ++neur) {
        
        std::string neurstr = Form("Layer_1_%i",neur);

        if (map_weight_h1_param_layer.at(paramstr).at(neurstr)) {
          content = fXML->GetNodeContent(map_weight_h1_param_layer.at(paramstr).at(neurstr));
          double d = stod(content);
          fWeight_h1(param,neur) = d;
        }
      }
    }
  }
  
  for (int neur=0; neur<fNNeuron; ++neur) {
    std::string neurstr = Form("Layer_1_%i",neur);
    content = fXML->GetNodeContent(map_weight_out.at(neurstr));
    double d = stod(content);
    fWeight_out(neur,0)=d;
  }
  
  for (int neur=0; neur<fNNeuron; ++neur) {
    std::string neurstr = Form("Layer_1_%i",neur);
    content = fXML->GetNodeContent(map_bias_b1.at(neurstr));
    double d = stod(content);
    fBias_b1(neur,0) = d;
  }

  for (auto const& it: map_bias_out) {
    content = fXML->GetNodeContent(it.second);
    double d = stod(content);
    fBias_out(0,0) = d;
  }

};


void ClusterEnergyEstimator::LoadChildAndAssertExistance(const XMLNodePointer_t node,
                                                         std::map<std::string,XMLNodePointer_t>& map_name_node,
                                                         const std::vector<std::string>& requiredchild,
                                                         const int nchild) const {
  
  std::cout << "Checking that node "  << fXML->GetNodeName(node) << " contains the child :" << std::endl;
  for (auto const& it:requiredchild)
    std::cout << " - " << it << std::endl;
  std::cout << nchild << " times." << std::endl;

  XMLNodePointer_t child = fXML->GetChild(node);
  while (child != 0) {
    map_name_node[fXML->GetNodeName(child)] = child;
    child = fXML->GetNext(child);
  }
  for (int i=0; i<nchild; ++i) {
    std::string addon = Form("_%i",i);
    if (nchild==1) addon="";
    for (auto const& it:requiredchild)
      if (map_name_node.find(it+addon) == map_name_node.end()) {
        std::cout << "XML File is corrupted!!" << std::endl;
        std::cout << "Node " << fXML->GetNodeName(node) << " doesn't contain the child " << it+addon << std::endl;
        std::cout << "Exiting..." << std::endl;
        exit(1);
      }
  }
  std::cout << "DONE Checking node "  << fXML->GetNodeName(node) << std::endl;
    
};



