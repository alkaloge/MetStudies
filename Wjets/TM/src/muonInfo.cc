#include "Wjets/TM/interface/muonInfo.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "Wjets/TM/interface/rhoInfo.h"

muonInfo::muonInfo(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  if(debug) std::cout<<"in muon constructor"<<std::endl;
  muonLabel_     = iConfig.getUntrackedParameter<edm::InputTag> ("muonLabel_");
  Muon_4Momentum                      = new TClonesArray("TLorentzVector");
  if(debug) std::cout<<"in muon constructor: calling SetBrances()"<<std::endl;
  SetBranches();
}

muonInfo::~muonInfo(){
  delete tree_;
  delete Muon_4Momentum;
}

void muonInfo::Fill(const edm::Event& iEvent, math::XYZPoint& pv, reco::Vertex& vtx, float& rhoLepton){
  Clear();
  if(debug_)    std::cout<<"getting muon info"<<std::endl;

  edm::Handle<edm::View<pat::Muon> > muonHandle;
  iEvent.getByLabel(muonLabel_,muonHandle);
  const edm::View<pat::Muon> & muons = *muonHandle;

  if(not iEvent.getByLabel(muonLabel_,muonHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "<<muonLabel_<<std::endl; 
    exit(0);
  }  

  edm::View<pat::Muon>::const_iterator muon;
  for(muon = muons.begin(); muon!=muons.end(); muon++){
    idPreselection = 0.;
    TLorentzVector p4(muon->px(),muon->py(),muon->pz(),muon->energy());
    new( (*Muon_4Momentum)[Muon_n]) TLorentzVector(p4);
    Muon_isLooseMuon.push_back((bool)muon->isLooseMuon());
    Muon_isMediumMuon.push_back((bool)muon->isMediumMuon());

    if(muon->innerTrack().isNonnull()){
      dxy = muon->innerTrack()->dxy(pv);
      dz = muon->innerTrack()->dz(pv);
    }


    float neutralHadIsoMu = muon->pfIsolationR04().sumNeutralHadronEt; 
    float photonIsoMu = muon->pfIsolationR04().sumPhotonEt;
    float chargedHadIsoMu = muon->pfIsolationR04().sumChargedHadronPt;
    float puIsoMu = muon->pfIsolationR04().sumPUPt;

    double neutralIsoMuN = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
    double neutralIsoMu = max(double(0),neutralIsoMuN); 
    float absIsoMu = chargedHadIsoMu + neutralIsoMu;
    float relIsoMu = absIsoMu/muon->pt();



    //if(muon->pt() > 10 && abs(muon->eta()) < 2.4 && abs(dxy) < 0.05 && abs(dz) < 0.1 && relIsoMu < 0.4 && muon->isLooseMuon()) idPreselection = 1.;
    if(muon->pt() > 10 && fabs(muon->eta()) < 2.4 && fabs(dxy) < 0.05 && fabs(dz) < 0.1 && relIsoMu < 0.4 && muon->isLooseMuon()) idPreselection = 1.;
    Muon_idPreselection.push_back(idPreselection);
    Muon_n++;
  }  
  if(debug_)    std::cout<<"got muon info"<<std::endl;
}

void muonInfo::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  AddBranch(&Muon_n  ,"Muon_n");
  AddBranch(&Muon_4Momentum ,"Muon_4Momentum");
  AddBranch(&Muon_isLooseMuon         ,"Muon_isLooseMuon");
  AddBranch(&Muon_isMediumMuon         ,"Muon_isMediumMuon");
  AddBranch(&Muon_idPreselection      ,"Muon_idPreselection");
  AddBranch(&relIsoMu      ,"relIsoMu");
  AddBranch(&dxy      ,"Muon_dxy");
  AddBranch(&dz      ,"Muon_dz");
  if(debug_)    std::cout<<"set branches"<<std::endl;
}

void muonInfo::Clear(){
  if(debug_)std::cout<<"clearing Muon info"<<std::endl;
  Muon_n = 0; 
  idPreselection = -99.;
  Muon_4Momentum->Clear();
  Muon_isLooseMuon.clear();
  Muon_isMediumMuon.clear();
  Muon_idPreselection.clear();
  if(debug_) std::cout<<"cleared"<<std::endl;
}
