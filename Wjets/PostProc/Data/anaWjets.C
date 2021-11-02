#define anaWjets_cxx
#include "anaWjets.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include "LepGood.h"
#include "XYMETCorrection.h"

std::pair<double,double> upara_uperp(double metx, double mety, double qx, double qy, double qt){
	double uX = - metx - qx;
	double uY = - mety - qy;
	double upara = (uX*qx + uY*qy)/qt;
	double uperp = (uX*qy - uY*qx)/qt;
	return {upara, uperp};
}

int main(int argc, char *argv[])
{

  if(argc > 1)
    {
      anaWjets t(argv[1], argv[2]);
      t.Loop();
    }
  return 0;
}

using namespace std;

void anaWjets::Loop()
{
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();
        Long64_t nbytes = 0, nb = 0;
        nentries=10000;
        string isomu1("HLT_IsoMu24_v"), isomu2("HLT_IsoMu27_v");

        size_t found1, found2;
        double upara, uperp, uX, uY, metX_corr, metY_corr;
        int prescale=9999; 
	int bin1, bin2, nJet40;
        vector<bool> clean_jet;
        vector<LepGood> veto_collection, LepGood_muon, LepGood_electron; LepGood loose_muon, loose_electron;
        bool tr24, tr27,  pass;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if (jentry%1000 == 0 ) cout<<" netry"<<jentry<<" "<<int(nentries)<<endl;

        upara = 0; uperp = 0; uX = 0; uY = 0; metX_corr = 0.0; metY_corr = 0.0;
        clean_jet.clear();
        
        LepGood_muon.clear(); LepGood_electron.clear(); veto_collection.clear(); nJet40 = 0; nvtx = 0;
        tr24 = false; tr27 = false; pass = false;
	jet_pt = -9999; jet_eta = -9999; fail = false;
	//bin1 = std::max(1, std::min(g2->GetNbinsX(), g2->GetXaxis()->FindBin((int)trueInteractions)));
//Trigger
	for(unsigned int t = 0; t < all_ifTriggerpassed->size(); t++){
                found1 = (*all_triggernames)[t].find(isomu1);
                found2 = (*all_triggernames)[t].find(isomu2);
                if(found1 != string::npos){ tr24 = (*all_ifTriggerpassed)[t]; }
                if(found2 != string::npos){ tr27 = (*all_ifTriggerpassed)[t]; }
        }
        if (tr24 || tr27) pass = true;
          
//Choose tight photons with Summer17 v2 id: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_94X_and_later
//Event selection
//Choose photon in barrel, add met filter requirements
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#2018_data
        if( (*filter_goodVertices)[0] && (*filter_globalsupertighthalo2016)[0] && (*filter_hbher2t)[0] && (*filter_hbheiso)[0] && (*filter_ecaltp)[0] && (*filter_ecalBadCalib)[0] && (*filter_badPFMuon)[0]){
                for(int i = 0; i < Vertex_n; i++){
                        TVector3* vtx = (TVector3*) (*Vertex_position)[i];
                        if(!(*Vertex_isFake)[i] && (*Vertex_ndof)[i] > 4 && abs(vtx->z()) <= 24. && abs((*Vertex_d0)[i]) <= 2.) nvtx++;
                }

                int good_muon=-1;
                for(int i = 0; i < Muon_n; i++){
                        TLorentzVector* mu = (TLorentzVector*) (*Muon_4Momentum)[i];
                        if((*Muon_idPreselection)[i] > 0.5 && mu->Pt()>29 ){
                                loose_muon.momentum = mu;
                                LepGood_muon.push_back(loose_muon);
                                good_muon=i;
                        }
                }

                if (LepGood_muon.size() != 1 ) continue;
                

                for(int i = 0; i < Electron_n; i++){
                        TLorentzVector* ele = (TLorentzVector*) (*Electron_4Momentum)[i];
                        int fail_ele = 0;
                        for (unsigned int k = 0; k < LepGood_muon.size(); k++){
                                if( (*LepGood_muon[k].momentum).DeltaR(*ele) < 0.2 ) 
                            fail_ele++;
                        }
                        if((*Electron_idPreselection)[i] > 0.5 && fail_ele > 0){
                                //loose_muon.momentum = ele;
                                //LepGood_coll.push_back(loose_electron);
                                pass = false;
                        }
                }
                if (pass==false) continue;

                for(int i = 0; i < PFJet_n; i++){
			TLorentzVector* j = (TLorentzVector*) (*PFJet_4Momentum)[i];
			if((j->Pt() > 40) && (abs(j->Eta()) <= 2.4) && ((*PFJet_NHF)[i] < 0.90) && ((*PFJet_NEMF)[i] < 0.90) && ((*PFJet_NumConst)[i] > 1) 
                                && ((*PFJet_CHF)[i] > 0) && ((*PFJet_CHM)[i] > 0)){

                        if ((*PFJet_DeepCSVb)[i] >=0.4184) pass = false;

                        for (unsigned int k = 0; k < LepGood_muon.size(); k++){
                                if( (*LepGood_muon[k].momentum).DeltaR(*j) > 0.5 ) 
             			clean_jet.push_back((bool) true);
			}
			}
			else clean_jet.push_back((bool) false); 
                }
                if (pass==false) continue;

//Select events using trigger

                for(int i = 0; i < PFJet_n; i++){
                        if(clean_jet[i]) nJet40++;
                }
                nJet = nJet40;

                TLorentzVector* q = (TLorentzVector*) (*Muon_4Momentum)[good_muon];

//Addtional quality criteria for photon, lepton veto
                //if(pass && LepGood_coll.size() == 0 && q->Pt() > 50 && abs(q->Eta()) < 1.44 && (*Photon_r9)[tight_photon[0]] > 0.9 && (*Photon_r9)[tight_photon[0]] < 1.0 && !(*Photon_hasPixelSeed)[tight_photon[0]] && (*Photon_passElectronVeto)[tight_photon[0]])
                        if (pass==true && LepGood_muon.size()>0){
                        //int etabin = std::max(1, std::min(h2->GetNbinsX(), h2->GetXaxis()->FindBin((q->Eta())));
                        //int ptbin = std::max(1, std::min(h2->GetNbinsY(), h2->GetYaxis()->FindBin(q->Pt())));
                        //gammaSF = h2->GetBinContent(etabin,ptbin);
                        //gammaSF = 1;
                        uparallel = upara_uperp((*PFMetPx)[0], (*PFMetPy)[0], q->Px(), q->Py(), q->Pt()).first;
			uperpendicular = upara_uperp((*PFMetPx)[0], (*PFMetPy)[0], q->Px(), q->Py(), q->Pt()).second;
                        lepton_pt = q->Pt();
                        lepton_eta = q->Eta();
                        lepton_phi = q->Phi();
                        lepton_energy = q->E();
                        met = (*PFMetPt)[0];
			met_phi = (*PFMetPhi)[0];
 
//std::pair<double,double> METXYCorr_Met_MetPhi(double uncormet, double uncormet_phi, int runnb, int year, bool isMC, int npv)
			corrmet = METXYCorr_Met_MetPhi(met, met_phi, RunNumber, 2018, false, Vertex_n).first;
			corrmet_phi = METXYCorr_Met_MetPhi(met, met_phi, RunNumber, 2018, false, Vertex_n).second;
			metX_corr = METXYCorr_Met_MetPhi(met, met_phi, RunNumber, 2018, false, Vertex_n).first;
			metY_corr = METXYCorr_Met_MetPhi(met, met_phi, RunNumber, 2018, false, Vertex_n).second;
                        corrupara = upara_uperp(metX_corr, metY_corr, q->Px(), q->Py(), q->Pt()).first;
                        corruperp = upara_uperp(metX_corr, metY_corr, q->Px(), q->Py(), q->Pt()).second;//

                        //wgt = prescale;
                        /*
			bin2 = std::max(1, std::min(g4->GetNbinsX(), g4->GetXaxis()->FindBin(nvtx)));
			nvsf = g4->GetBinContent(bin2);
			wgtA = gammaSF*(g6->GetBinContent(bin2));
                        wgtB = gammaSF*(g8->GetBinContent(bin2));
                        wgtC = gammaSF*(g10->GetBinContent(bin2));
                        wgtD = gammaSF*(g12->GetBinContent(bin2));
                        */
                        for(int i = 0; i < PFJet_n; i++){
                                if(clean_jet[i]){
                                	TLorentzVector* recjet = (TLorentzVector*) (*PFJet_4Momentum)[i];
                                	jet_pt = recjet->Pt();
                                	jet_eta = recjet->Eta();
					jet_phi = recjet->Phi();
                                	break;
                                }
                        }
			//IDSF = gammaSF;
                        mettree->Fill();
                }
	} 
   }
}

void anaWjets::BookHistos(const char* file2){
        fileName = new TFile(file2, "RECREATE");
        fileName->cd();
        char name[100];

        mettree = new TTree("mettree", "selected events");
        mettree->Branch("uparallel",&uparallel,"uparallel/F");
        mettree->Branch("uperpendicular",&uperpendicular,"uperpendicular/F");
        mettree->Branch("lepton_pt",&lepton_pt,"lepton_pt/F");
	mettree->Branch("lepton_eta",&lepton_eta,"lepton_eta/F");
        mettree->Branch("lepton_phi",&lepton_phi,"lepton_phi/F");
        mettree->Branch("lepton_energy",&lepton_energy,"lepton_energy/F");
        mettree->Branch("met",&met,"met/F");
        mettree->Branch("met_phi",&met_phi,"met_phi/F");
        mettree->Branch("corrmet",&corrmet,"corrmet/F");
        mettree->Branch("corrmet_phi",&corrmet_phi,"corrmet_phi/F");
        mettree->Branch("corrupara",&corrupara,"corrupara/F");
        mettree->Branch("corruperp",&corruperp,"corruperp/F");
	mettree->Branch("nvsf",&nvsf,"nvsf/F");
	mettree->Branch("nvtx",&nvtx,"nvtx/I");
	mettree->Branch("nJet", &nJet, "nJet/I");
        mettree->Branch("jet_pt",&jet_pt,"jet_pt/F");
        mettree->Branch("jet_eta",&jet_eta,"jet_eta/F");
        mettree->Branch("jet_phi",&jet_phi,"jet_phi/F");
        mettree->Branch("fail",&fail,"fail/B");
}
