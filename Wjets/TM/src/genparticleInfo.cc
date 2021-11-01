#include "Wjets/TM/interface/genparticleInfo.h"

genparticleInfo::genparticleInfo(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  if(debug) std::cout<<"in genparticle constructor"<<std::endl;
  genParticleLabel_               = iConfig.getUntrackedParameter<edm::InputTag> ("genParticleLabel_");
  Gen_Muon_4Momentum            = new TClonesArray("TLorentzVector");

  if(debug) std::cout<<"in genparticle constructor: calling SetBrances()"<<std::endl;
  SetBranches();
}

genparticleInfo::~genparticleInfo(){
  delete tree_;
  delete Gen_Muon_4Momentum            ;
}

void genparticleInfo::Fill(const edm::Event& iEvent){
  Clear();
  if(debug_)    std::cout<<"getting genparticle info"<<std::endl;
  
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(genParticleLabel_, genParticles);
  if(not iEvent.getByLabel(genParticleLabel_, genParticles )){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "<<genParticleLabel_<<std::endl;
    exit(0);
  }

  reco::GenParticleCollection genColl(*(genParticles.product()));
  std::sort(genColl.begin(),genColl.end(),PtGreater());
  for (reco::GenParticleCollection::const_iterator genparticle = genParticles->begin(); genparticle != genParticles->end(); genparticle++){
      if (abs(genparticle->pdgId()) == 13 && genparticle->status() == 1){ 
	  new( (*Gen_Muon_4Momentum)[Gen_Muon_n]) TLorentzVector((float)genparticle->px(),
								     (float)genparticle->py(),
								     (float)genparticle->pz(),
								     (float)genparticle->energy());
	  Gen_Muon_isPrompt.push_back(genparticle->isPromptFinalState());
	  Gen_Muon_status.push_back(genparticle->status());
	  Gen_Muon_n++;
	}
    }
 
  if(debug_) cout<<"got genparticle info"<<endl;  
}


void genparticleInfo::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  AddBranch(&Gen_Muon_n                    ,"Gen_Muon_n");
  AddBranch(&Gen_Muon_isPrompt             ,"Gen_Muon_isPrompt");
  AddBranch(&Gen_Muon_status            ,"Gen_Muon_status");
  AddBranch(&Gen_Muon_4Momentum            ,"Gen_Muon_4Momentum");
 
  if(debug_)    std::cout<<"set branches"<<std::endl;
}

void genparticleInfo::Clear(){
  if(debug_)std::cout<<"clearing Muon info"<<std::endl;

  Gen_Muon_n = 0;
  Gen_Muon_isPrompt.clear();
  Gen_Muon_4Momentum->Clear();            

  if(debug_) std::cout<<"cleared"<<std::endl;
}

