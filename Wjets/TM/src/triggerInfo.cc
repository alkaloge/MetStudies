#include "GammaJets/TM/interface/triggerInfo.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

triggerInfo::triggerInfo(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  if(debug) std::cout<<"in triggerInfo constructor"<<std::endl;
  HLTriggerResults_ = iConfig.getUntrackedParameter<edm::InputTag>("HLTriggerResults_");
  ntriggers = 0;
  if(debug) std::cout<<"in trigger constructor: calling SetBrances()"<<std::endl;
  SetBranches();
}

triggerInfo::~triggerInfo(){
  delete tree_;
}


unsigned int triggerInfo::AddTriggerObjects(const edm::Event& iEvent, 
					    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> TriggerObjectCollectionToken,
					    const edm::TriggerResults & triggerResults) {

  // trigger objects
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(TriggerObjectCollectionToken, triggerObjects);
  
  assert(triggerObjects.isValid());

for (unsigned int iTO=0; iTO<triggerObjects->size(); ++iTO) {
    if (trigobject_count==M_trigobjectmaxcount) {
      cerr << "number of trigger objects > M_trigobjectmaxcount. They are missing." << endl; 
      break;
    }

   //(*triggerObjects)[iTO].unpackNamesAndLabels( iEvent, triggerResults  );
    //    std::vector<std::string> filterLabels = (*triggerObjects)[iTO].pathNames();
    std::vector<std::string> filterLabels = (*triggerObjects)[iTO].filterLabels();
    bool matchFound = false;
    //    std::cout << "Object " << iTO << " filters : " <<  filterLabels.size() << std::endl;
    //    for (unsigned int i=0; i < filterLabels.size(); ++i) {
    //      TString FilterLabel(filterLabels.at(i));
    //      if (FilterLabel.Contains("hltL1sMu18erTau24erIorMu20erTau24er"))
    //	std::cout << "    " << filterLabels.at(i) << std::endl;
    //    }
    std::vector<bool> passedFilters; passedFilters.clear();
    //    std::vector<std::string> matchedFilters; matchedFilters.clear();
    for (unsigned int n=0; n < run_hltfilters.size(); ++n) {
      TString HltFilter(run_hltfilters.at(n));
      bool thisMatched = false;
      for (unsigned int i=0; i < filterLabels.size(); ++i) {
	TString FilterName(filterLabels.at(i));
	if (HltFilter==FilterName) {
	  matchFound = true;
	  thisMatched = true;
	  //	  matchedFilters.push_back(filterLabels.at(i));
	  break;
	}
      }
      passedFilters.push_back(thisMatched);
    } 
    if (matchFound) {
      //      std::cout << "   trigger object " << iTO 
      //		<< "   pt = " << (*triggerObjects)[iTO].pt() 
      //		<< "   eta = " << (*triggerObjects)[iTO].eta() 
      //		<< "   phi = " << (*triggerObjects)[iTO].phi() << std::endl;
      //      for (unsigned int ifilter=0; ifilter<matchedFilters.size(); ++ifilter)
      //	std::cout << "    " << matchedFilters[ifilter] << std::endl;
      for (unsigned int n=0; n < M_hltfiltersmax; ++n) {
	if (n<passedFilters.size())
	  trigobject_filters[trigobject_count][n] = passedFilters.at(n);
	else
	  trigobject_filters[trigobject_count][n] = false;
      }
      trigobject_px[trigobject_count] = (*triggerObjects)[iTO].px();
      trigobject_py[trigobject_count] = (*triggerObjects)[iTO].py();
      trigobject_pz[trigobject_count] = (*triggerObjects)[iTO].pz();
      trigobject_pt[trigobject_count] = (*triggerObjects)[iTO].pt();
      trigobject_eta[trigobject_count] = (*triggerObjects)[iTO].eta();
      trigobject_phi[trigobject_count] = (*triggerObjects)[iTO].phi();
      trigobject_isMuon[trigobject_count] = (*triggerObjects)[iTO].hasTriggerObjectType(trigger::TriggerMuon);
      trigobject_isElectron[trigobject_count] = (*triggerObjects)[iTO].hasTriggerObjectType(trigger::TriggerElectron);
      trigobject_isTau[trigobject_count] = (*triggerObjects)[iTO].hasTriggerObjectType(trigger::TriggerTau);
      trigobject_isJet[trigobject_count] = (*triggerObjects)[iTO].hasTriggerObjectType(trigger::TriggerJet);
      trigobject_isMET[trigobject_count] = (*triggerObjects)[iTO].hasTriggerObjectType(trigger::TriggerMET);
      trigobject_count++;
    }
  }



return trigobject_count;
}


void triggerInfo::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup, std::vector<std::string>& all_triggers, HLTConfigProvider& hltConfig_, HLTPrescaleProvider& hltPrescale_, std::string& hltlabel_, const size_t& MaxN =200 ){

  if(debug_)    std::cout<<"getting HL Trigger info"<<std::endl;

  Clear();
  edm::Handle<TriggerResults> HLTR;
  iEvent.getByLabel(HLTriggerResults_,HLTR);
  if(debug_) std::cout<<"prescale set: "<<hltPrescale_.prescaleSet( iEvent, iSetup)<<std::endl;

   if (HLTR.isValid())
     {
       const edm::TriggerNames &triggerNames_ = iEvent.triggerNames(*HLTR);
       
       vector<int> idx;
      Int_t hsize = Int_t(HLTR->size());


      ntriggers = all_triggers.size();
      for(int i = 0; i< ntriggers;i++){

  	all_triggerprescales.push_back(0);
  	all_ifTriggerpassed.push_back(0);
  	idx.push_back(triggerNames_.triggerIndex(all_triggers[i]));
 	
  	if(idx[i] < hsize){
          all_triggernames.push_back(all_triggers[i]);
  	  all_ifTriggerpassed[i]=HLTR->accept(idx[i]);
          all_triggerprescales[i] = hltPrescale_.prescaleValue( iEvent, iSetup, all_triggers[i]);
 	}	
 	
  	if(debug_){
  	  std::cout<<"prescale for "<<all_triggers[i]<<" is: "<< all_triggerprescales[i]<<std::endl;
  	  std::cout<<"if triggger passed for "<<all_triggers[i]<<" : "<<all_ifTriggerpassed[i]<<std::endl;
  	}
	
      }//loop over trigger
      
    }//HLT is valid collection
  
  if(debug_)    std::cout<<"got trigger info"<<std::endl;
}

void triggerInfo::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  AddBranch(&all_triggernames,"all_triggernames");
  AddBranch(&all_triggerprescales,"all_triggerprescales");
  AddBranch(&all_ifTriggerpassed ,"all_ifTriggerpassed" );
    if(debug_)    std::cout<<"set branches"<<std::endl;
}

void triggerInfo::Clear(){
  if(debug_)std::cout<<"clearing trigger info"<<std::endl;
  all_triggernames.clear();
  all_triggerprescales.clear();
  all_ifTriggerpassed.clear();
  if(debug_) std::cout<<"cleared"<<std::endl;
}
