#include "Wjets/TM/interface/genparticleInfo.h"


using namespace reco;
using namespace std;

genparticleInfo::genparticleInfo(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector cc):
//genparticleInfo::genparticleInfo(std::string name, TTree* tree, bool debug, const pset& iConfig):
  baseTree(name,tree,debug),
  //genEventInfoToken_(cc.consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfo")))
   GenToken_( cc.consumes<GenEventInfoProduct>(edm::InputTag("generator")))
{
  if(debug) std::cout<<"in genparticle constructor"<<std::endl;
  genParticleLabel_               = iConfig.getUntrackedParameter<edm::InputTag> ("genParticleLabel_");
  Gen_Muon_4Momentum            = new TClonesArray("TLorentzVector");
  Gen_Electron_4Momentum            = new TClonesArray("TLorentzVector");

  if(debug) std::cout<<"in genparticle constructor: calling SetBrances()"<<std::endl;
  SetBranches();
}

genparticleInfo::~genparticleInfo(){
  delete tree_;
  delete Gen_Muon_4Momentum            ;
  delete Gen_Electron_4Momentum            ;
}

void genparticleInfo::Fill(const edm::Event& iEvent,TH1D*h1) {
  Clear();
  if(debug_)    std::cout<<"getting genparticle info"<<std::endl;
  
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(genParticleLabel_, genParticles);

  if(not iEvent.getByLabel(genParticleLabel_, genParticles )){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "<<genParticleLabel_<<std::endl;
    exit(0);
  }


  genWeight=1;
  edm::Handle<GenEventInfoProduct> GenEventInfo;
  iEvent.getByToken( GenToken_, GenEventInfo);
  if(GenEventInfo.isValid())
  {
   genWeight = GenEventInfo->weight();
   //h1->Fill(0,genweight);
   h1->Fill(0.,genWeight);
	}



  reco::GenParticleCollection genColl(*(genParticles.product()));
  std::sort(genColl.begin(),genColl.end(),PtGreater());
  for (reco::GenParticleCollection::const_iterator genparticle = genParticles->begin(); genparticle != genParticles->end(); genparticle++){
      if (abs(genparticle->pdgId()) == 13 && genparticle->status() == 1){ 
      //if (abs(genparticle->pdgId()) == 13 ){ 
	  new( (*Gen_Muon_4Momentum)[Gen_Muon_n]) TLorentzVector((float)genparticle->px(),
								     (float)genparticle->py(),
								     (float)genparticle->pz(),
								     (float)genparticle->energy());
	  Gen_Muon_isPrompt.push_back(genparticle->isPromptFinalState());
	  Gen_Muon_status.push_back(genparticle->status());
	  Gen_Muon_n++;
   //cout << "Event weight for muon  : " << genparticle->isHardProcess() << endl;
	}
      if (abs(genparticle->pdgId()) == 11 && genparticle->status() == 1){ 
      //if (abs(genparticle->pdgId()) == 13 ){ 
	  new( (*Gen_Electron_4Momentum)[Gen_Electron_n]) TLorentzVector((float)genparticle->px(),
								     (float)genparticle->py(),
								     (float)genparticle->pz(),
								     (float)genparticle->energy());
	  Gen_Electron_isPrompt.push_back(genparticle->isPromptFinalState());
	  Gen_Electron_status.push_back(genparticle->status());
	  Gen_Electron_n++;
   //cout << "Event weight for muon  : " << genparticle->isHardProcess() << endl;
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
  AddBranch(&Gen_Electron_n                    ,"Gen_Electron_n");
  AddBranch(&Gen_Electron_isPrompt             ,"Gen_Electron_isPrompt");
  AddBranch(&Gen_Electron_status            ,"Gen_Electron_status");
  AddBranch(&Gen_Electron_4Momentum            ,"Gen_Electron_4Momentum");
  AddBranch(&genWeight            ,"genWeight");
 
  if(debug_)    std::cout<<"set branches"<<std::endl;
}

void genparticleInfo::Clear(){
  if(debug_)std::cout<<"clearing Muon info"<<std::endl;

  Gen_Muon_n = 0;
  Gen_Muon_isPrompt.clear();
  Gen_Muon_4Momentum->Clear();            
  Gen_Muon_status.clear();            
  Gen_Electron_n = 0;
  Gen_Electron_isPrompt.clear();
  Gen_Electron_4Momentum->Clear();            
  Gen_Electron_status.clear();            
  genWeight=0;

  if(debug_) std::cout<<"cleared"<<std::endl;
}

