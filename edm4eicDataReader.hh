#ifndef edm4eicDataReader_hh
#define edm4eicDataReader_hh 

#include <ROOT/RDataFrame.hxx>
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"
#include "TH1D.h"
#include "TH2D.h"

#include <fstream> 
#include <iostream>
#include <memory>
#include <functional>
#include <thread>
#include <string>
#include <vector>
#include <unordered_map>
#include <numbers>

unsigned int nThreads = 5;

std::unordered_map<std::string, std::vector<std::shared_ptr<TH1D>>> hist1DMaps;
std::unordered_map<std::string, std::vector<std::shared_ptr<TH2D>>> hist2DMaps;
std::vector<double> jetPtBins = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 120};
std::vector<double> jetEtaBins = {-2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1
                            , 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5}; 
std::vector<double> jetPhiBins = {-3.2, -3.0, -2.8, -2.6, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2
                            , 0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2};
std::vector<double> dRBins = {};

std::vector<std::vector<std::string>> fileLists(nThreads);
std::vector<std::ofstream> logs(nThreads);

float R = 1.0;
float conPtMin = 0.2;
float conPtMax = 100;

//bool doReadJets = false;
bool doClusterJets = true;
bool doBijectiveMatching = true;

#endif

#ifdef edm4eicDataReader_cpp
template<typename T>
double getDeltaR(T& jet1, T& jet2){
    double dEta = std::abs(jet1.eta() - jet2.eta());
    double dPhi = std::abs(jet1.phi() - jet2.phi());
    if(dPhi > std::numbers::pi)dPhi = 2.*std::numbers::pi - dPhi;
    return std::sqrt(dEta*dEta + dPhi*dPhi);
}
#endif