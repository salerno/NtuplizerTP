#ifndef Ntuplizer_H
#define Ntuplizer_H

// CMSSW
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <DataFormats/HepMCCandidate/interface/GenParticleFwd.h>
#include "DataFormats/Math/interface/LorentzVector.h"
#include <FWCore/Framework/interface/ESHandle.h>

// ROOT
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

// C++
#include<memory>
#include<vector>
#include <map>

//#include "Analyzer/Ntuplizer/interface/LinkDef.h"

using namespace std;

class Ntuplizer : public edm::EDAnalyzer {
   public:
      explicit Ntuplizer(const edm::ParameterSet&);
      ~Ntuplizer();

      typedef math::XYZTLorentzVector LorentzVector ;
      
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      void Init();
      void FillEvent(const edm::Event&, const edm::EventSetup&);
      void FillElectrons(const edm::Event&, const edm::EventSetup&);
      void FillVertices(const edm::Event&, const edm::EventSetup&);
      void FillTruth(const edm::Event&, const edm::EventSetup&);
      void FillMET (const edm::Event& iEvent, const edm::EventSetup& iSetup);
      void FillTracks (const edm::Event& iEvent, const edm::EventSetup& iSetup);
      
      //void setMomentum(TLorentzVector & myvector, const LorentzVector & mom);
      void setMomentum(TLorentzVector & myvector, const LorentzVector & mom) ;


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      //inputTag
      edm::InputTag EleTag_;
      edm::InputTag VerticesTag_;
      edm::InputTag TracksTag_;
      bool isMC_;	
      edm::InputTag PileupSrc_;	
      
	
      //tree
      TTree *_mytree;
      TLorentzVector myvector ;  

      //global variables
      int _nEvent, _nRun, _nLumi;
      //pile-up
      int _PU_N;
      //vertices
      int _vtx_N;
      //rho
      float _rho;
      
      // MET
      float _met_pf_et, _met_pf_px, _met_pf_py, _met_pf_phi, _met_pf_set, _met_pf_sig; 

      //electrons
      int ele_N;
      TClonesArray * m_electrons ;

      //std vector
      std::vector<int> ele_echarge ;
      //energy matches
      std::vector<float> ele_he ; 
      std::vector<float> ele_hebc ; 
      std::vector<float> ele_eseedpout ; 
      std::vector<float> ele_ep ; 
      std::vector<float> ele_eseedp ; 
      std::vector<float> ele_eelepout ;       
      //delta 
      std::vector<float> ele_deltaetaseed ; 
      std::vector<float> ele_deltaetaele ; 
      std::vector<float> ele_deltaphiseed ; 
      std::vector<float> ele_deltaphiele ; 
      std::vector<float> ele_deltaetain ; 
      std::vector<float> ele_deltaphiin ;
      //cluster shapes
      std::vector<float> ele_sigmaietaieta ; 
      std::vector<float> ele_sigmaetaeta ; 
      std::vector<float> ele_sigmaiphiiphi; 
      std::vector<float> ele_e15 ; 
      std::vector<float> ele_e25max ; 
      std::vector<float> ele_e55 ; 
      std::vector<float> ele_r9 ;
      //momentum
      std::vector<float> ele_pin_mode ; 
      std::vector<float> ele_pout_mode ; 
      std::vector<float> ele_pTin_mode ; 
      std::vector<float> ele_pTout_mode ; 
      //brem
      std::vector<float> ele_fbrem ;
      std::vector<float> ele_SCfbrem ;
      std::vector<float> ele_pfSCfbrem ; 
      std::vector<float> ele_trackfbrem ;
      std::vector<int> ele_nbrem ;
      //position
      std::vector<int> ele_isbarrel ; 
      std::vector<int> ele_isendcap ; 
      std::vector<int> ele_isEBetaGap ; 
      std::vector<int> ele_isEBphiGap ; 
      std::vector<int> ele_isEEdeeGap ; 
      std::vector<int> ele_isEEringGap ;
      std::vector<int> ele_isecalDriven ; 
      std::vector<int> ele_istrackerDriven ;
      std::vector<int> ele_eClass ;
      //distance vertex
      //std::vector<float> ele_dxyB ; 
      std::vector<float> ele_dxy ; 
      //std::vector<float> ele_dzB ; 
      std::vector<float> ele_dz ; 
      //std::vector<float> ele_dszB ; 
      std::vector<float> ele_dsz ;              
      //conversion
      std::vector<int> ele_valid_hits ; 
      std::vector<int> ele_lost_hits ; 
      std::vector<int> ele_chi2_hits ; 
      std::vector<float> ele_conv_dcot ;
      std::vector<float> ele_conv_dist ;
      std::vector<float> ele_conv_radius ;
      std::vector<int> ele_expected_inner_hits ;
      //mva 
      std::vector<float>  ele_mva ;
      //isolation
      std::vector<float> ele_pfChargedHadIso ; 
      std::vector<float> ele_pfNeutralHadIso ; 
      std::vector<float> ele_pfPhotonIso ; 
      //supercluster energies
      std::vector<float> ele_sclE ; 
      std::vector<float> ele_sclEt ; 
      std::vector<float> ele_sclEta ; 
      std::vector<float> ele_sclPhi ; 
      std::vector<float> ele_sclRawE ;
      std::vector<int>   ele_sclNclus ;
      //supercluster variables
      std::vector<float>  ele_sclphiwidth ; 
      std::vector<float>  ele_scletawidth ;
      //errors
      std::vector<float> ele_ecalErr ; 
      std::vector<float> ele_trackErr ; 
      std::vector<float> ele_combErr ; 
      std::vector<float> ele_PFcombErr;
      //kf       
      std::vector<float>  ele_kfchi2 ;
      std::vector<int>  ele_kfhits ;
      
      //pf variables
      std::vector<int> ele_pf_number;
      std::vector<int> ele_pf_id;
      std::vector<float> ele_pf_eta;
      std::vector<float> ele_pf_phi; 
      std::vector<float> ele_pf_pt; 
      std::vector<float> ele_pf_dz; 
      std::vector<float> ele_pf_dxy;
      std::vector<float> ele_pf_vx; 
      std::vector<float> ele_pf_vy; 
      std::vector<float> ele_pf_vz;
      std::vector<float> ele_pf_mva_nog; 
      std::vector<float> ele_pf_mva_epi;

      /*
      std::vector<int>> ele_pf_id;
      std::vector<float>> ele_pf_eta;
      std::vector<float>> ele_pf_phi; 
      std::vector<float>> ele_pf_pt; 
      std::vector<float>> ele_pf_dz; 
      std::vector<float>> ele_pf_dxy;
      std::vector<float>> ele_pf_vx; 
      std::vector<float>> ele_pf_vy; 
      std::vector<float>> ele_pf_vz;
      std::vector<float>> ele_pf_mva_nog; 
      std::vector<float>> ele_pf_mva_epi;
      */
      
      std::vector<int>  ele_seed_subDet2 ;
      std::vector<float> ele_seed_dRz2 ;
      std::vector<float> ele_seed_dPhi2 ;
      std::vector<float> ele_seed_dRz2Pos ;
      std::vector<float> ele_seed_dPhi2Pos ;
      std::vector<int> ele_seed_subDet1 ; 
      std::vector<float> ele_seed_dRz1;
      std::vector<float> ele_seed_dPhi1;
      std::vector<float> ele_seed_dRz1Pos;
      std::vector<float> ele_seed_dPhi1Pos;
      
	TClonesArray * _m_MC_gen_V;
	TClonesArray * _m_MC_gen_Higgs;
	TClonesArray * _m_MC_gen_leptons;
	TClonesArray * _m_MC_gen_leptons_status1;
	TClonesArray * _m_MC_gen_leptons_status2;
	float _MC_gen_V_pdgid[10];
	float _MC_gen_leptons_pdgid[10];
	float _MC_gen_leptons_status1_pdgid[10];
	float _MC_gen_leptons_status2_pdgid[10];
	int _MC_flavor[2];

};
#endif
