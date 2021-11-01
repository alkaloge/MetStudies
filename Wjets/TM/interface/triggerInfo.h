#ifndef __TRIGGER_INFO_H_
#define __TRIGGER_INFO_H_


#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "Wjets/TM/interface/utils.h"
#include "Wjets/TM/interface/baseTree.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
using namespace std;
using namespace edm;
using namespace reco;
#define M_trigobjectmaxcount 1000
#define M_hltfiltersmax 200

class triggerInfo : public baseTree{

 public:
  triggerInfo(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg);
  ~triggerInfo();
  void Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup, 
	    std::vector<std::string>& all_triggers, HLTConfigProvider& hltConfig_, HLTPrescaleProvider& hltPrescale_,
	    std::string& hltlabel_, const size_t& MaxN );
  unsigned int AddTriggerObjects(const edm::Event& iEvent, edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> TriggerObjectCollectionToken, const edm::TriggerResults & trigRes);

  void SetBranches();
  void Clear();

 private:
  triggerInfo(){};
  edm::InputTag HLTriggerResults_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> TriggerObjectCollectionToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> myTriggerObjectCollectionToken_;
  vector<string> run_hltfilters;
  edm::Handle<trigger::TriggerEvent> HLTriggerEvent;
  UInt_t trigobject_count;
  Float_t trigobject_px[M_trigobjectmaxcount];
  Float_t trigobject_py[M_trigobjectmaxcount];
  Float_t trigobject_pz[M_trigobjectmaxcount];
  Float_t trigobject_pt[M_trigobjectmaxcount];
  Float_t trigobject_eta[M_trigobjectmaxcount];
  Float_t trigobject_phi[M_trigobjectmaxcount];
  Bool_t trigobject_isMuon[M_trigobjectmaxcount];
  Bool_t trigobject_isElectron[M_trigobjectmaxcount];
  Bool_t trigobject_isTau[M_trigobjectmaxcount];
  Bool_t trigobject_isJet[M_trigobjectmaxcount];
  Bool_t trigobject_isMET[M_trigobjectmaxcount];
  Bool_t  trigobject_filters[M_trigobjectmaxcount][M_hltfiltersmax];

  int ntriggers;
  std::vector<int>  all_triggerprescales;
  std::vector<bool> all_ifTriggerpassed;
  std::vector<std::string> all_triggernames; 
};

#endif

