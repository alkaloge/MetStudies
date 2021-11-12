#ifndef __GENPARTICLE_INFO_H_
#define __GENPARTICLE_INFO_H_
#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include <set>
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "Wjets/TM/interface/utils.h"
#include "Wjets/TM/interface/baseTree.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TMath.h"
#include "TTree.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


using namespace std;
using namespace reco;

class genparticleInfo : public baseTree
{

 public:
  genparticleInfo(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector cc);
 // genparticleInfo(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg);
  ~genparticleInfo();
  void Fill(const edm::Event& iEvent, TH1D*h1);
  void SetBranches();
  void Clear();
  
 private:
  genparticleInfo(){};
  edm::InputTag genParticleLabel_;
  //edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_; //desy
  edm::EDGetTokenT<GenEventInfoProduct > GenToken_;

  int Gen_Muon_n;
  TClonesArray *Gen_Muon_4Momentum;
  std::vector<bool> Gen_Muon_isPrompt; 
  std::vector<int> Gen_Muon_status; 
  int Gen_Electron_n;
  TClonesArray *Gen_Electron_4Momentum;
  std::vector<bool> Gen_Electron_isPrompt; 
  std::vector<int> Gen_Electron_status; 
  float genWeight; 
};

#endif

