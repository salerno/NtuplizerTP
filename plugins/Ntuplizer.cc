// MY include
#include "Ntuplizer.h"

// C++ include files
#include <memory>

// CMSSW include
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/PhotonCore.h"
#include "DataFormats/EgammaCandidates/interface/PhotonCoreFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexSorter.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//
// class declaration
//

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//using namespace std;
using namespace reco;
using namespace edm;

// =============================================================================================
// constructor
Ntuplizer::Ntuplizer(const edm::ParameterSet& iConfig) :
// ==============================================================================================
EleTag_ (iConfig.getParameter<edm::InputTag> ("EleTag")),
PhoTag_ (iConfig.getParameter<edm::InputTag> ("PhoTag")),
VerticesTag_(iConfig.getParameter<edm::InputTag> ("VerticesTag")),
TracksTag_(iConfig.getParameter<edm::InputTag> ("TracksTag")),
isMC_ (iConfig.getParameter<bool>("isMC")),
PileupSrc_ ("addPileupInfo")
{
}

// =============================================================================================
// destructor
Ntuplizer::~Ntuplizer()
// =============================================================================================
{
  delete m_electrons;
  delete m_photons;
  
  if(isMC_ ) {
    delete _m_MC_gen_V;
    delete _m_MC_gen_leptons;
    delete _m_MC_gen_Higgs;
    delete _m_MC_gen_leptons_status1;
    delete _m_MC_gen_leptons_status2;
    delete _m_MC_gen_photonsHiggs;
    delete _m_MC_gen_photonsPrompt;
  } // if MC
  
  
  
}

// =============================================================================================
// ------------ method called once each job just before starting event loop  ------------
void Ntuplizer::beginJob()
//=============================================================================================
{
  // Book histograms
  edm::Service<TFileService> fs;
  _mytree  = fs->make <TTree>("simpleRootTree","simpleRootTree"); 
  
  // Global
  _mytree->Branch("nEvent",&_nEvent,"nEvent/I");
  _mytree->Branch("nRun",&_nRun,"nRun/I");
  // Pile UP
  _mytree->Branch("PU_N",&_PU_N,"PU_N/I");
  // Vertices
  _mytree->Branch("vtx_N",&_vtx_N,"vtx_N/I");
  // rho
  _mytree->Branch("rho",&_rho,"rho/F");
  // Electrons
  _mytree->Branch("ele_N",&ele_N,"ele_N/I");
  m_electrons = new TClonesArray ("TLorentzVector");
  _mytree->Branch ("electrons", "TClonesArray", &m_electrons, 256000,0);
  
      _mytree->Branch("ele_echarge", &ele_echarge );
      //energy matches
      _mytree->Branch("ele_he", &ele_he ); 
      _mytree->Branch("ele_hebc", &ele_hebc ); 
      _mytree->Branch("ele_eseedpout", &ele_eseedpout ); 
      _mytree->Branch("ele_ep", &ele_ep ); 
      _mytree->Branch("ele_eseedp", &ele_eseedp ); 
      _mytree->Branch("ele_eelepout", &ele_eelepout );
      //momentum
      _mytree->Branch("ele_pin_mode", &ele_pin_mode ); 
      _mytree->Branch("ele_pout_mode", &ele_pout_mode ); 
      _mytree->Branch("ele_pTin_mode", &ele_pTin_mode ); 
      _mytree->Branch("ele_pTout_mode", &ele_pTout_mode ); 
      //delta 
      _mytree->Branch("ele_deltaetaseed", &ele_deltaetaseed ); 
      _mytree->Branch("ele_deltaetaele", &ele_deltaetaele ); 
      _mytree->Branch("ele_deltaphiseed", &ele_deltaphiseed ); 
      _mytree->Branch("ele_deltaphiele", &ele_deltaphiele ); 
      _mytree->Branch("ele_deltaetain", &ele_deltaetain ); 
      _mytree->Branch("ele_deltaphiin", &ele_deltaphiin );
      //cluster shapes
      _mytree->Branch("ele_sigmaietaieta", &ele_sigmaietaieta ); 
      _mytree->Branch("ele_sigmaetaeta", &ele_sigmaetaeta ); 
      _mytree->Branch("ele_sigmaiphiiphi", &ele_sigmaiphiiphi ); 
      _mytree->Branch("ele_e15", &ele_e15 ); 
      _mytree->Branch("ele_e25max", &ele_e25max ); 
      _mytree->Branch("ele_e55", &ele_e55 ); 
      _mytree->Branch("ele_r9", &ele_r9 );
      //brem
      _mytree->Branch("ele_fbrem", &ele_fbrem );
      _mytree->Branch("ele_SCfbrem", &ele_SCfbrem );
      _mytree->Branch("ele_pfSCfbrem", &ele_pfSCfbrem ); 
      _mytree->Branch("ele_trackfbrem", &ele_trackfbrem );
      _mytree->Branch("ele_nbrem", &ele_nbrem );
       //position
      _mytree->Branch("ele_isbarrel", &ele_isbarrel ); 
      _mytree->Branch("ele_isendcap", &ele_isendcap ); 
      _mytree->Branch("ele_isEBetaGap", &ele_isEBetaGap ); 
      _mytree->Branch("ele_isEBphiGap", &ele_isEBphiGap ); 
      _mytree->Branch("ele_isEEdeeGap", &ele_isEEdeeGap ); 
      _mytree->Branch("ele_isEEringGap", &ele_isEEringGap );
      _mytree->Branch("ele_isecalDriven", &ele_isecalDriven ); 
      _mytree->Branch("ele_istrackerDriven", &ele_istrackerDriven );
      _mytree->Branch("ele_eClass", &ele_eClass );
      //errors
      _mytree->Branch("ele_ecalErr", &ele_ecalErr ); 
      _mytree->Branch("ele_trackErr", &ele_trackErr ); 
      _mytree->Branch("ele_combErr", &ele_combErr ); 
      _mytree->Branch("ele_PFcombErr", &ele_PFcombErr ); 
      //mva 
      _mytree->Branch("ele_mva", &ele_mva );
      //tracking Variables
      _mytree->Branch("ele_valid_hits", &ele_valid_hits ); 
      _mytree->Branch("ele_lost_hits", &ele_lost_hits ); 
      _mytree->Branch("ele_chi2_hits", &ele_chi2_hits ); 
      _mytree->Branch("ele_kfchi2", &ele_kfchi2 );
      _mytree->Branch("ele_kfhits", &ele_kfhits );
      //position wrt beam or 000 
      //_mytree->Branch("ele_dxyB", &ele_dxyB); 
      _mytree->Branch("ele_dxy", &ele_dxy ); 
      //_mytree->Branch("ele_dzB", &ele_dzB ); 
      _mytree->Branch("ele_dz", &ele_dz ); 
      //_mytree->Branch("ele_dszB", &ele_dszB ); 
      _mytree->Branch("ele_dsz", &ele_dsz );              
      //conversion      
      _mytree->Branch("ele_conv_dcot", &ele_conv_dcot );
      _mytree->Branch("ele_conv_dist", &ele_conv_dist );
      _mytree->Branch("ele_conv_radius", &ele_conv_radius );
      _mytree->Branch("ele_expected_inner_hits", &ele_expected_inner_hits );
      //isolation
      _mytree->Branch("ele_pfChargedHadIso", &ele_pfChargedHadIso ); 
      _mytree->Branch("ele_pfNeutralHadIso", &ele_pfNeutralHadIso ); 
      _mytree->Branch("ele_pfPhotonIso", &ele_pfPhotonIso ); 
       //supercluster 
      _mytree->Branch("ele_sclE", &ele_sclE ); 
      _mytree->Branch("ele_sclEt", &ele_sclEt ); 
      _mytree->Branch("ele_sclEta", &ele_sclEta ); 
      _mytree->Branch("ele_sclPhi", &ele_sclPhi); 
      _mytree->Branch("ele_sclRawE", &ele_sclRawE );
      _mytree->Branch("ele_sclNclus", &ele_sclNclus );
      _mytree->Branch("ele_sclphiwidth", &ele_sclphiwidth ); 
      _mytree->Branch("ele_scletawidth", &ele_scletawidth );

//seeds
  _mytree->Branch("ele_seed_subDet2", &ele_seed_subDet2);
  _mytree->Branch("ele_seed_dRz2",&ele_seed_dRz2);
  _mytree->Branch("ele_seed_dPhi2",&ele_seed_dPhi2);
  _mytree->Branch("ele_seed_dRz2Pos", &ele_seed_dRz2Pos);
  _mytree->Branch("ele_seed_dPhi2Pos",&ele_seed_dPhi2Pos);
  _mytree->Branch("ele_seed_subDet1",&ele_seed_subDet1); 
  _mytree->Branch("ele_seed_dRz1",&ele_seed_dRz1);
  _mytree->Branch("ele_seed_dPhi1",&ele_seed_dPhi1);
  _mytree->Branch("ele_seed_dRz1Pos",&ele_seed_dRz1Pos);
  _mytree->Branch("ele_seed_dPhi1Pos", &ele_seed_dPhi1Pos);

   _mytree->Branch("ele_pf_number", &ele_pf_number);
   _mytree->Branch("ele_pf_id", &ele_pf_id);
   _mytree->Branch("ele_pf_phi", &ele_pf_phi);
   _mytree->Branch("ele_pf_eta", &ele_pf_eta);
   _mytree->Branch("ele_pf_pt", &ele_pf_pt);
   _mytree->Branch("ele_pf_dz", &ele_pf_dz);
   _mytree->Branch("ele_pf_dxy", &ele_pf_dxy);
   _mytree->Branch("ele_pf_vx", &ele_pf_vx);
   _mytree->Branch("ele_pf_vy", &ele_pf_vy); 
   _mytree->Branch("ele_pf_vz", &ele_pf_vz);
   _mytree->Branch("ele_pf_mva_nog", &ele_pf_mva_nog); 
   _mytree->Branch("ele_pf_mva_epi", &ele_pf_mva_epi);
   
// Photons
  _mytree->Branch("pho_N",&pho_N,"pho_N/I");
  m_photons = new TClonesArray ("TLorentzVector");
  _mytree->Branch ("photons", "TClonesArray", &m_photons, 256000,0);

	_mytree->Branch("pho_pt", &pho_pt );
	_mytree->Branch("pho_phi", &pho_phi ); 
	_mytree->Branch("pho_eta", &pho_eta );  
	//position
	_mytree->Branch("pho_isEB", &pho_isEB );  
	_mytree->Branch("pho_isEE", &pho_isEE );
	_mytree->Branch("pho_isEBEEGap", &pho_isEBEEGap );    
	_mytree->Branch("pho_isEBGap", &pho_isEBGap );  
	_mytree->Branch("pho_isEBEtaGap", &pho_isEBEtaGap );  
	_mytree->Branch("pho_isEBPhiGap", &pho_isEBPhiGap );
	_mytree->Branch("pho_isEEGap", &pho_isEEGap );
	_mytree->Branch("pho_isEEDeeGap", &pho_isEEDeeGap );
	_mytree->Branch("pho_isEERingGap", &pho_isEERingGap );
	//Shower Shape
	_mytree->Branch("pho_sigmaietaieta", &pho_sigmaietaieta );  
	_mytree->Branch("pho_sigmaetaeta", &pho_sigmaetaeta );  
	_mytree->Branch("pho_e15", &pho_e15 );
	_mytree->Branch("pho_e25", &pho_e25 );
	_mytree->Branch("pho_e33", &pho_e33 );
	_mytree->Branch("pho_e55", &pho_e55 );
	_mytree->Branch("pho_r9", &pho_r9 );
	_mytree->Branch("pho_maxEnergyXtal", &pho_maxEnergyXtal );
	//SuperCluster
	_mytree->Branch("pho_sclEtaWidth", &pho_sclEtaWidth );
	_mytree->Branch("pho_sclPhiWidth", &pho_sclPhiWidth );
	_mytree->Branch("pho_sclRawE", &pho_sclRawE );
	//HoE
	_mytree->Branch("pho_EcalRechitsSumEt", &pho_EcalRechitsSumEt );  
	_mytree->Branch("pho_HcalTowerSumEt", &pho_HcalTowerSumEt );
	_mytree->Branch("pho_hoe", &pho_hoe );
	//PF isolation
        _mytree->Branch("pho_pf_number", &pho_pf_number);
        _mytree->Branch("pho_pf_id", &pho_pf_id);
        _mytree->Branch("pho_pf_phi", &pho_pf_phi);
        _mytree->Branch("pho_pf_eta", &pho_pf_eta);
        _mytree->Branch("pho_pf_pt", &pho_pf_pt);
        _mytree->Branch("pho_pf_dz", &pho_pf_dz);
        _mytree->Branch("pho_pf_dxy", &pho_pf_dxy);
        _mytree->Branch("pho_pf_vz", &pho_pf_vz);	
        _mytree->Branch("pho_pf_vy", &pho_pf_vy);			
        _mytree->Branch("pho_pf_vx", &pho_pf_vx);					
        _mytree->Branch("pho_pf_mva_nog", &pho_pf_mva_nog); 
        _mytree->Branch("pho_pf_mva_epi", &pho_pf_mva_epi);
	//Photon Object Definition
        _mytree->Branch("pho_isPFlowPhoton",&pho_isPFlowPhoton); 
        _mytree->Branch("pho_isStandardPhoton",&pho_isStandardPhoton);
	//Electron Veto
	_mytree->Branch("pho_hasPixelSeed", &pho_hasPixelSeed );


  // PFMET
  _mytree->Branch("met_pf_et",&_met_pf_et,"met_pf_et/F");
  _mytree->Branch("met_pf_px",&_met_pf_px,"met_pf_px/F");
  _mytree->Branch("met_pf_py",&_met_pf_py,"met_pf_py/F");
  _mytree->Branch("met_pf_phi",&_met_pf_phi,"met_pf_phi/F");
  _mytree->Branch("met_pf_set",&_met_pf_set,"met_pf_set/F");
  _mytree->Branch("met_pf_sig",&_met_pf_sig,"met_pf_sig/F");

  // Truth 
  _m_MC_gen_V = new TClonesArray ("TLorentzVector");
  _mytree->Branch ("MC_gen_V", "TClonesArray", &_m_MC_gen_V, 256000,0);
  _mytree->Branch ("MC_gen_V_pdgid",&_MC_gen_V_pdgid, "MC_gen_V_pdgid[10]/F");
  //
  _m_MC_gen_Higgs = new TClonesArray ("TLorentzVector");
  _mytree->Branch ("MC_gen_Higgs", "TClonesArray", &_m_MC_gen_Higgs, 256000,0);
  //
  _m_MC_gen_leptons         = new TClonesArray ("TLorentzVector");
  _m_MC_gen_leptons_status1 = new TClonesArray ("TLorentzVector");
  _m_MC_gen_leptons_status2 = new TClonesArray ("TLorentzVector");
  _mytree->Branch ("MC_gen_leptons", "TClonesArray", &_m_MC_gen_leptons, 256000,0);
  _mytree->Branch ("MC_gen_leptons_status1", "TClonesArray", &_m_MC_gen_leptons_status1, 256000,0);
  _mytree->Branch ("MC_gen_leptons_status2", "TClonesArray", &_m_MC_gen_leptons_status2, 256000,0);
  _mytree->Branch ("MC_gen_leptons_pdgid",&_MC_gen_leptons_pdgid, "MC_gen_leptons_pdgid[10]/F");
  _mytree->Branch ("MC_gen_leptons_status1_pdgid",&_MC_gen_leptons_status1_pdgid, "MC_gen_leptons_status1_pdgid[10]/F");
  _mytree->Branch ("MC_gen_leptons_status2_pdgid",&_MC_gen_leptons_status2_pdgid, "MC_gen_leptons_status2_pdgid[10]/F");
  _mytree->Branch ("MC_flavor",&_MC_flavor,"MC_flavor[2]/I");
  //
  _m_MC_gen_photonsHiggs = new TClonesArray ("TLorentzVector");
  _mytree->Branch ("MC_gen_photonsHiggs", "TClonesArray", &_m_MC_gen_photonsHiggs, 256000,0);
  _m_MC_gen_photonsPrompt = new TClonesArray ("TLorentzVector");
  _mytree->Branch ("MC_gen_photonsPrompt", "TClonesArray", &_m_MC_gen_photonsPrompt, 256000,0);

}


//
// member functions
//
// =============================================================================================
// ------------ method called for each event  ------------
void Ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// =============================================================================================
{
   using namespace edm;

   Init();

   FillEvent(iEvent, iSetup);
   FillVertices(iEvent, iSetup);
   FillTracks(iEvent, iSetup);
   m_electrons -> Clear();
   FillElectrons(iEvent, iSetup);
   m_photons -> Clear();
   FillPhotons(iEvent, iSetup);
   FillMET (iEvent, iSetup);

   if(isMC_ ) {
     _m_MC_gen_V->Clear();
     _m_MC_gen_Higgs->Clear();
     _m_MC_gen_leptons->Clear();
     _m_MC_gen_leptons_status1->Clear();
     _m_MC_gen_leptons_status2->Clear();
     _m_MC_gen_photonsHiggs->Clear();
     _m_MC_gen_photonsPrompt->Clear();
     FillTruth(iEvent, iSetup);
   }

   _mytree->Fill();

// #ifdef THIS_IS_AN_EVENT_EXAMPLE
//    Handle<ExampleData> pIn;
//    iEvent.getByLabel("example",pIn);
// #endif
   
// #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
//    ESHandle<SetupData> pSetup;
//    iSetup.get<SetupRecord>().get(pSetup);
// #endif
}


// =============================================================================================
void Ntuplizer::FillTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup)
//=============================================================================================
{
   Handle<vector<reco::Track> >  recoGeneralTracksCollection;
   iEvent.getByLabel(TracksTag_, recoGeneralTracksCollection);
   /*
   for(reco::Track::const_iterator itTrack = recoGeneralTracksCollection->begin(); 
   itTrack != recoGeneralTracksCollection->end(); 
       ++itTrack) {
     std::cout<<"looping on tracks"<<std::endl;
   }
   */
} // end of FillEvent


// =============================================================================================
void Ntuplizer::FillEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup)
//=============================================================================================
{

   Handle<vector<PileupSummaryInfo> > PupInfo;
   iEvent.getByLabel(PileupSrc_, PupInfo);
	for (vector<PileupSummaryInfo>::const_iterator cand = PupInfo->begin();cand != PupInfo->end(); ++ cand) {
			_PU_N = cand->getPU_NumInteractions();
	} // loop on Pile up

  _nEvent = iEvent.id().event();
  _nRun   = iEvent.id().run();
  _nLumi  = iEvent.luminosityBlock();
  
  
  InputTag theRhoTag = InputTag("kt6PFJets","rho","RECO");
  edm::Handle<double> rhoHandle;
  iEvent.getByLabel(theRhoTag, rhoHandle);
  double rho = *rhoHandle;
  _rho = rho ;  
 
} // end of FillEvent



// =============================================================================================
void Ntuplizer::FillVertices(const edm::Event& iEvent, const edm::EventSetup& iSetup)
//=============================================================================================
{
  
  Handle<vector<reco::Vertex> >  recoPrimaryVertexCollection;
  iEvent.getByLabel(VerticesTag_, recoPrimaryVertexCollection);
  _vtx_N = recoPrimaryVertexCollection->size();
  
} // end of FillVertices

// =============================================================================================
void Ntuplizer::FillElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup)
//=============================================================================================
{
  
  edm::Handle<reco::GsfElectronCollection> electronsCol;
  iEvent.getByLabel(EleTag_, electronsCol);
  
  InputTag  vertexLabel(string("offlinePrimaryVertices"));
  Handle<reco::VertexCollection> thePrimaryVertexColl;
  iEvent.getByLabel(VerticesTag_ ,thePrimaryVertexColl);
  
  edm::Handle< reco::PFCandidateCollection > pfHandle;
  iEvent.getByLabel("particleFlow", pfHandle);
  const reco::PFCandidateCollection* pfCandColl = pfHandle.product();
  
  TClonesArray & electrons = *m_electrons;
  int counter = 0;
  ele_N = electronsCol->size();

  for (reco::GsfElectronCollection::const_iterator ielectrons=electronsCol->begin(); ielectrons!=electronsCol->end();++ielectrons) {
    if(counter>49) continue;

    setMomentum(myvector, ielectrons->p4());
    new (electrons[counter]) TLorentzVector(myvector);

    //here loop over the PF candidates around the electrons
    int numberOfPFAroundEle = 0;
    for(unsigned i=0; i<pfCandColl->size(); i++) {
      const reco::PFCandidate& pfc = (*pfCandColl)[i];
      float dR = deltaR(ielectrons->eta(), ielectrons->phi(), pfc.momentum().Eta(), pfc.momentum().Phi());
      int pfID = pfc.particleId();
      if ( (dR < 0.4) && (pfID == 1 || pfID == 4 || pfID == 5) ) 
      {
       ++numberOfPFAroundEle;
       float pf_pt =  pfc.pt();
       float pf_eta = pfc.momentum().Eta();
       float pf_phi = pfc.momentum().Phi();
	   float pf_mva_nog = pfc.mva_nothing_gamma();
       float pf_mva_epi = pfc.mva_e_pi();
       float dz_pf =  -999.;
       float dxy_pf =   -999.;
       float vx_pf =  -999.;
       float vy_pf =  -999.;
       float vz_pf =  -999.;
       if (pfID == 1) { //|| pfID == 2 || pfID == 3) { // charged hadrons || electrons || muons 
       dz_pf =  pfc.trackRef()->dz();
       //float dz_ele =  ielectrons->gsfTrack()->dz();
       dxy_pf =  pfc.trackRef()->dxy();
       //float dxy_ele =  ielectrons->gsfTrack()->dxy();
       vx_pf = pfc.vx();
       vy_pf = pfc.vy();
       vz_pf = pfc.vz();
       }
           ele_pf_id.push_back(pfID);
           ele_pf_eta.push_back(pf_eta); 
           ele_pf_phi.push_back(pf_phi); 
           ele_pf_pt.push_back(pf_pt); 
           ele_pf_dz.push_back(dz_pf); 
           ele_pf_dxy.push_back(dxy_pf); 
           ele_pf_vx.push_back(vx_pf);
           ele_pf_vy.push_back(vy_pf);
           ele_pf_vz.push_back(vz_pf);
           ele_pf_mva_nog.push_back(pf_mva_nog); 
           ele_pf_mva_epi.push_back(pf_mva_epi);
     }
    }
    ele_pf_number.push_back(numberOfPFAroundEle);

    //
    ele_echarge.push_back(ielectrons->charge()); 
    //
    ele_he.push_back(ielectrons->hcalOverEcal()); 
    ele_hebc.push_back(ielectrons-> hcalOverEcalBc());
    // TrackCluster Matching
    ele_eseedpout.push_back(ielectrons->eSeedClusterOverPout());
    ele_ep.push_back(ielectrons->eSuperClusterOverP());        
    ele_eseedp.push_back(ielectrons->eSeedClusterOverP());         
    ele_eelepout.push_back(ielectrons->eEleClusterOverPout());       
    //
    ele_pin_mode.push_back(ielectrons->trackMomentumAtVtx().R()); 
    ele_pout_mode.push_back(ielectrons->trackMomentumOut().R()); 
    ele_pTin_mode.push_back(ielectrons->trackMomentumAtVtx().Rho()); 
    ele_pTout_mode.push_back(ielectrons->trackMomentumOut().Rho()); 
    //
    ele_deltaetaseed.push_back(ielectrons->deltaEtaSeedClusterTrackAtCalo()); 
    ele_deltaphiseed.push_back(ielectrons->deltaPhiSeedClusterTrackAtCalo());  
    ele_deltaetaele.push_back(ielectrons->deltaEtaEleClusterTrackAtCalo());  
    ele_deltaphiele.push_back(ielectrons->deltaPhiEleClusterTrackAtCalo()); 
    ele_deltaetain.push_back(ielectrons->deltaEtaSuperClusterTrackAtVtx());
    ele_deltaphiin.push_back(ielectrons->deltaPhiSuperClusterTrackAtVtx());   
    // Shower Shape
    ele_sigmaietaieta.push_back((ielectrons->showerShape()).sigmaIetaIeta); 
    ele_sigmaetaeta.push_back((ielectrons->showerShape()).sigmaEtaEta); 
    ele_sigmaiphiiphi.push_back((ielectrons->showerShape()).sigmaIphiIphi);
    ele_e15.push_back((ielectrons->showerShape()).e1x5);
    ele_e25max.push_back((ielectrons->showerShape()).e2x5Max);
    ele_e55.push_back((ielectrons->showerShape()).e5x5);
    ele_r9.push_back((ielectrons->showerShape()).r9); 
    // fbrem
    ele_fbrem.push_back(ielectrons->fbrem());
    ele_trackfbrem.push_back(ielectrons->trackFbrem());
    ele_SCfbrem.push_back(ielectrons->superClusterFbrem());
    ele_pfSCfbrem.push_back(ielectrons->pfSuperClusterFbrem()); // Should be identical to the previous one...
    ele_nbrem.push_back(ielectrons->numberOfBrems());
    //position
    ele_eClass.push_back(ielectrons->classification());
    if (ielectrons->isEB()) ele_isbarrel.push_back(1); 
    else  ele_isbarrel.push_back(0);
    if (ielectrons->isEE()) ele_isendcap.push_back(1); 
    else  ele_isendcap.push_back(0);
    if (ielectrons->isEBEtaGap()) ele_isEBetaGap.push_back(1);  
    else ele_isEBetaGap.push_back(0);  
    if (ielectrons->isEBPhiGap()) ele_isEBphiGap.push_back(1); 
    else  ele_isEBphiGap.push_back(0); 
    if (ielectrons->isEEDeeGap()) ele_isEEdeeGap.push_back(1);  
    else ele_isEEdeeGap.push_back(0);
    if (ielectrons->isEERingGap()) ele_isEEringGap.push_back(1);
    else ele_isEEringGap.push_back(0);
    if (ielectrons->ecalDrivenSeed()) ele_isecalDriven.push_back(1);
    else ele_isecalDriven.push_back(0);
    if (ielectrons->trackerDrivenSeed()) ele_istrackerDriven.push_back(1);
    else ele_istrackerDriven.push_back(0);
    // E/P combination
    ele_ecalErr.push_back(ielectrons->ecalEnergyError());
    ele_trackErr.push_back(ielectrons->trackMomentumError());
    ele_combErr.push_back(ielectrons->p4Error(GsfElectron::P4_COMBINATION));
    ele_PFcombErr.push_back(ielectrons->p4Error(GsfElectron::P4_PFLOW_COMBINATION));
    // mva
    ele_mva.push_back(ielectrons->mva());
    // Tracking Variables
    ele_lost_hits.push_back(ielectrons->gsfTrack()->lost()); 
    ele_valid_hits.push_back(ielectrons->gsfTrack()->found()); 
    ele_chi2_hits.push_back(ielectrons->gsfTrack()->normalizedChi2());
    bool validKF=false;
    reco::TrackRef myTrackRef = ielectrons->closestCtfTrackRef();
    validKF = myTrackRef.isNonnull();
    ele_kfchi2.push_back(validKF ? myTrackRef->normalizedChi2() : 0); 
    ele_kfhits.push_back(validKF ? myTrackRef->hitPattern().trackerLayersWithMeasurement() : -1.); 
    //position wrt beam or 000 
    //ele_dxyB.push_back(ielectrons->gsfTrack()->dxy(bs.position()));
    ele_dxy.push_back(ielectrons->gsfTrack()->dxy());
    //ele_dzB.push_back(ielectrons->gsfTrack()->dz(bs.position()));
    ele_dz.push_back(ielectrons->gsfTrack()->dz());
    //ele_dszB.push_back(ielectrons->gsfTrack()->dsz(bs.position()));
    ele_dsz.push_back(ielectrons->gsfTrack()->dsz());
    //conversion
    ele_conv_dcot.push_back(ielectrons->convDist()); 
    ele_conv_dist.push_back(ielectrons->convDcot()); 
    ele_conv_radius.push_back(ielectrons->convRadius());
    ele_expected_inner_hits.push_back(ielectrons->gsfTrack()->trackerExpectedHitsInner().numberOfHits());
    // PF Isolation
    ele_pfChargedHadIso.push_back((ielectrons->pfIsolationVariables()).chargedHadronIso);
    ele_pfNeutralHadIso.push_back((ielectrons->pfIsolationVariables()).neutralHadronIso);
    ele_pfPhotonIso.push_back((ielectrons->pfIsolationVariables()).photonIso);
    //supercluster
    reco::SuperClusterRef sclRef = ielectrons->superCluster();
    double R  = TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y() +sclRef->z()*sclRef->z());
    double Rt = TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y());
    ele_sclRawE.push_back(sclRef->rawEnergy());
    ele_sclE.push_back(sclRef->energy());
    ele_sclEt.push_back(sclRef->energy()*(Rt/R));
    ele_sclEta.push_back(sclRef->eta());
    ele_sclPhi.push_back(sclRef->phi());
    ele_sclNclus.push_back(sclRef->clustersSize());
    ele_sclphiwidth.push_back(sclRef->phiWidth());
    ele_scletawidth.push_back(sclRef->etaWidth());

    //reco::GsfTrackRef track = eleIter->gsfTrack();
    if(ielectrons->gsfTrack()->extra().isAvailable() && 
      ielectrons->gsfTrack()->extra()->seedRef().isAvailable()) {
      ElectronSeedRef seedRef =  ielectrons->gsfTrack()->extra()->seedRef().castTo<ElectronSeedRef>();
      //ElectronSeedRef seedRef =  ielectrons->gsfTrack()->seedRef().castTo<ElectronSeedRef>();

      ele_seed_subDet2.push_back(seedRef->subDet2())  ;
      ele_seed_dRz2.push_back(seedRef->dRz2())  ;
      ele_seed_dPhi2.push_back(seedRef->dPhi2())  ;
      ele_seed_dRz2Pos.push_back(seedRef->dRz2Pos())  ;
      ele_seed_dPhi2Pos.push_back(seedRef->dPhi2Pos())  ;
      ele_seed_subDet1.push_back(seedRef->subDet1())  ; 
      ele_seed_dRz1.push_back(seedRef->dRz1()) ;
      ele_seed_dPhi1.push_back(seedRef->dPhi1()) ;
      ele_seed_dRz1Pos.push_back(seedRef->dRz1Pos()) ;
      ele_seed_dPhi1Pos.push_back(seedRef->dPhi1Pos()) ;

    }

    ++counter;
  } // for loop on gsfelectrons

} // end of FillElectrons


// =============================================================================================
void Ntuplizer::FillPhotons(const edm::Event& iEvent, const edm::EventSetup& iSetup)
//=============================================================================================
{
  
  edm::Handle<reco::PhotonCollection> photonsCol;
  iEvent.getByLabel(PhoTag_, photonsCol);
  
  InputTag  vertexLabel(string("offlinePrimaryVertices"));
  Handle<reco::VertexCollection> thePrimaryVertexColl;
  iEvent.getByLabel(VerticesTag_ ,thePrimaryVertexColl);
  
  edm::Handle< reco::PFCandidateCollection > pfHandle;
  iEvent.getByLabel("particleFlow", pfHandle);
  const reco::PFCandidateCollection* pfCandColl = pfHandle.product();
  
  TClonesArray & photons = *m_photons;
  int counter = 0;
  pho_N = photonsCol->size();

  for (reco::PhotonCollection::const_iterator iphotons=photonsCol->begin(); iphotons!=photonsCol->end();++iphotons) {
    if(counter>49) continue;

    setMomentum(myvector, iphotons->p4());
    new (photons[counter]) TLorentzVector(myvector); 

    pho_eta.push_back(iphotons->eta());   
    pho_phi.push_back(iphotons->phi()); 
    pho_pt.push_back(iphotons->pt());
    //Shower Shape
    pho_sigmaietaieta.push_back(iphotons->sigmaIetaIeta());   
    pho_sigmaetaeta.push_back(iphotons->sigmaEtaEta());   
    pho_e15.push_back(iphotons->e1x5() );
    pho_e25.push_back(iphotons->e2x5() );
    pho_e33.push_back(iphotons->e3x3() );
    pho_e55.push_back(iphotons->e5x5() );
    pho_r9.push_back(iphotons->r9() );
    pho_maxEnergyXtal.push_back(iphotons->maxEnergyXtal());
    //SuperCluster
    pho_sclEtaWidth.push_back(iphotons->superCluster()->etaWidth());
    pho_sclPhiWidth.push_back(iphotons->superCluster()->phiWidth());
    pho_sclRawE.push_back(iphotons->superCluster()->rawEnergy()); 
    //Position     
    pho_isEE.push_back(iphotons->isEE());      
    pho_isEB.push_back(iphotons->isEB()); 
    pho_isEBEEGap.push_back(iphotons->isEBEEGap()); 
    pho_isEBGap.push_back(iphotons->isEBGap()); 
    pho_isEBEtaGap.push_back(iphotons->isEBEtaGap()); 
    pho_isEBPhiGap.push_back(iphotons->isEBPhiGap()); 
    pho_isEEGap.push_back(iphotons->isEEGap()); 
    pho_isEEDeeGap.push_back(iphotons->isEEDeeGap()); 
    pho_isEERingGap.push_back(iphotons->isEERingGap());   
    //HoE 
    pho_EcalRechitsSumEt.push_back(iphotons->ecalRecHitSumEtConeDR03());  
    pho_HcalTowerSumEt.push_back(iphotons->hcalTowerSumEtConeDR03() );
    pho_hoe.push_back(iphotons->hadronicOverEm());
    //Object Definition
    pho_isPFlowPhoton.push_back(iphotons->isPFlowPhoton());
    pho_isStandardPhoton.push_back(iphotons->isStandardPhoton());
    //Electron Veto        
    pho_hasPixelSeed.push_back(iphotons->hasPixelSeed() );

   //here loop over the PF candidates around the photons
    int numberOfPFAroundPho = 0;
    for(unsigned i=0; i<pfCandColl->size(); i++) {
      const reco::PFCandidate& pfc = (*pfCandColl)[i];
      float dR = deltaR(iphotons->eta(), iphotons->phi(), pfc.momentum().Eta(), pfc.momentum().Phi());
      int pfID = pfc.particleId();
      if ( (dR < 0.4) && (pfID == 1 || pfID == 2 || pfID == 3 || pfID==4 || pfID == 5) ) //charged Hadrons || electrons || muons || photons || neutral hadrons)
      {
       ++numberOfPFAroundPho;
       float pf_pt =  pfc.pt();
       float pf_eta = pfc.momentum().Eta();
       float pf_phi = pfc.momentum().Phi();
       float pf_mva_nog = pfc.mva_nothing_gamma();
       float pf_mva_epi = pfc.mva_e_pi();
       pho_pf_id.push_back(pfID);
       pho_pf_eta.push_back(pf_eta); 
       pho_pf_phi.push_back(pf_phi); 
       pho_pf_pt.push_back(pf_pt); 
       pho_pf_mva_nog.push_back(pf_mva_nog); 
       pho_pf_mva_epi.push_back(pf_mva_epi);
       float dz_pf = -999.;
       float dxy_pf = -999.;
       float vx_pf = -999.;
       float vy_pf = -999.;
       float vz_pf = -999.;
       if (pfID == 1) { 
	 dz_pf = pfc.trackRef()->dz();
	 dxy_pf = pfc.trackRef()->dxy();
	 vx_pf = pfc.vx();
	 vy_pf = pfc.vy();
	 vz_pf = pfc.vz();
       }
       pho_pf_dxy.push_back(dxy_pf);
       pho_pf_dz.push_back(dz_pf);
       pho_pf_vx.push_back(vx_pf);
       pho_pf_vy.push_back(vy_pf);
       pho_pf_vz.push_back(vz_pf);

     }
    }
    pho_pf_number.push_back(numberOfPFAroundPho);
    //

    ++counter;

  } // for loop on photons

} // end of FillPhotons


// ====================================================================================
void Ntuplizer::FillMET (const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
	
  edm::Handle< edm::View<reco::PFMET> > pfMEThandle;
  iEvent.getByLabel("pfMet", pfMEThandle);
  
  
  // PFMET
  _met_pf_et  = (pfMEThandle->front() ).et();
  _met_pf_px  = (pfMEThandle->front() ).px();
  _met_pf_py  = (pfMEThandle->front() ).py();
  _met_pf_phi = (pfMEThandle->front() ).phi();
  _met_pf_set = (pfMEThandle->front() ).sumEt();
  _met_pf_sig = (pfMEThandle->front() ).mEtSig();
  
} // end of Fill MET



// ====================================================================================
void Ntuplizer::FillTruth(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
  
  edm::Handle<View<Candidate> > genCandidatesCollection;
  iEvent.getByLabel("genParticles", genCandidatesCollection); //genParticlesPruned
  
  TClonesArray &MC_gen_V               = *_m_MC_gen_V;
  TClonesArray &MC_gen_Higgs           = *_m_MC_gen_Higgs;
  TClonesArray &MC_gen_leptons         = *_m_MC_gen_leptons;
  TClonesArray &MC_gen_leptons_status2 = *_m_MC_gen_leptons_status2;
  TClonesArray &MC_gen_leptons_status1 = *_m_MC_gen_leptons_status1;
  TClonesArray &MC_gen_photonsHiggs    = *_m_MC_gen_photonsHiggs;
  TClonesArray &MC_gen_photonsPrompt    = *_m_MC_gen_photonsPrompt;  

  
  int counter             = 0;
  int counter_higgs       = 0;
  int counter_daughters   = 0;
  int counter_lep_status2 = 0;
  int counter_lep_status1 = 0;
  int counter_photonsHiggs = 0;
  int counter_photonsPrompt  = 0;
  
  // ----------------------------
  //      Loop on particles
  // ----------------------------
  for( View<Candidate>::const_iterator p = genCandidatesCollection->begin();p != genCandidatesCollection->end(); ++ p ) {
    
    // %%%%%%%%%%%%%%%%%%
    // If Leptons from Z
    // %%%%%%%%%%%%%%%%%%
    if(fabs(p->pdgId())==11 || fabs(p->pdgId())==13 ||  fabs(p->pdgId())==15) {
      
      if(p->status()==1) {
	if(p->numberOfMothers()>0) { // Need a mother...
	  if(p->mother(0)->pdgId()== p->pdgId()) {
	    setMomentum(myvector, p->p4());
	    new (MC_gen_leptons_status1[counter_lep_status1]) TLorentzVector(myvector);
	    _MC_gen_leptons_status1_pdgid[counter_lep_status1] = p->pdgId();
	    counter_lep_status1++;
	  }// if pdgid
	} // if mother
      } // if status 1
      
      if(p->status()==3) {
	if(p->numberOfMothers()>0) { // Need a mother...
	  if(p->mother(0)->pdgId()==23) {  // If Mother is a Z 
	    if(p->numberOfDaughters()>0) { // Need a daughter...
	      // Status 2 Leptons
	      if(p->daughter(0)->pdgId()==p->pdgId() && p->daughter(0)->status()==2) { // if daughter is lepton & status 2
		setMomentum(myvector, p->daughter(0)->p4());
		new (MC_gen_leptons_status2[counter_lep_status2]) TLorentzVector(myvector);
		_MC_gen_leptons_status2_pdgid[counter_lep_status2] = p->daughter(0)->pdgId();
		counter_lep_status2++;
		  } // if Daughter Status = 2
	    } // Need a daughter
	  } // if MOther = Z
	} // if Nmother (status3) >0
      } // if hard scatter electron (status3)
    } // if leptons
    
    // %%%%%%%%%%%%%%%%%%
    //     If W or Z
    // %%%%%%%%%%%%%%%%%%
    if (p->pdgId() == 23 || fabs(p->pdgId())==24) {
      
      
      if(p->status()==3) {
	// Fill truth W,Z
	setMomentum (myvector,p->p4());
	new (MC_gen_V[counter]) TLorentzVector(myvector);
	_MC_gen_V_pdgid[counter] = p->pdgId();
	// Loop on daughters
	for(unsigned int i=0;i<p->numberOfDaughters();i++) {
	  bool islep = false;
	  if(fabs(p->daughter(i)->pdgId())==11) { _MC_flavor[counter] = 0; islep=true;} // electron
	  if(fabs(p->daughter(i)->pdgId())==13) { _MC_flavor[counter] = 1; islep=true;} // muon
	  if(fabs(p->daughter(i)->pdgId())==15) { _MC_flavor[counter] = 2; islep=true;} // taus
	  
	  if(islep) { // p->daughter(i)->status()==1) { ?!
	    setMomentum(myvector, p->daughter(i)->p4());
	    new (MC_gen_leptons[counter_daughters]) TLorentzVector(myvector);
	    _MC_gen_leptons_pdgid[counter_daughters] = p->daughter(i)->pdgId();
	    
	    counter_daughters++;
	  } // if is lepton
	} // for loop on daughters
	counter++;
      } // if status stable
    } // if W or Z


    // %%%%%%%%%%%%%%%%%%
    // If Higgs
    // %%%%%%%%%%%%%%%%%%
	if (p->pdgId() == 25 && p->status()==3) { //if Higgs
        	setMomentum (myvector,p->p4());
        	new (MC_gen_Higgs[counter_higgs]) TLorentzVector(myvector);
        	counter_higgs++;
        	// Loop on daughters (photons)
		for(unsigned int i=0;i<p->numberOfDaughters();i++) {
	  		if(fabs(p->daughter(i)->pdgId())==22) { 
	    			setMomentum(myvector, p->daughter(i)->p4());
	    			new (MC_gen_photonsHiggs[counter_photonsHiggs]) TLorentzVector(myvector);
	    			counter_photonsHiggs++;
	  		} 
		} 
        } 
    // %%%%%%%%%%%%%%%%%%
    // If prompt photons
    // %%%%%%%%%%%%%%%%%%
   	else if(p->numberOfMothers()>0) { // Need a mother...
   	if (p->pdgId() == 22 && p->mother(0)->pdgId()==22 && p->status()==1) { //if photon not produced by jet
        	setMomentum (myvector,p->p4());
      		new (MC_gen_photonsPrompt[counter_photonsPrompt]) TLorentzVector(myvector);
        	counter_photonsPrompt++;
	}
   	}

  } // for loop on particles
} // end of FillTruth

// ------------ method called once each job just after ending the event loop  ------------
void 
Ntuplizer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
Ntuplizer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
Ntuplizer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
Ntuplizer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
Ntuplizer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/


 // ====================================================================================================
void Ntuplizer::setMomentum(TLorentzVector &myvector, const LorentzVector & mom)
// ====================================================================================================
{
  
  myvector.SetPx (mom.Px());
  myvector.SetPy (mom.Py());
  myvector.SetPz (mom.Pz());
  myvector.SetE (mom.E());
  
}


 // ====================================================================================================
void Ntuplizer::Init()
// ====================================================================================================
{

  _PU_N = 0;
  
  _vtx_N = 0;

  ele_N = 0;
  
  _rho = 0.;
  
  _met_pf_et  = 0.;
  _met_pf_px  = 0.; 
  _met_pf_py  = 0.; 
  _met_pf_phi = 0.; 
  _met_pf_set = 0.; 
  _met_pf_sig = 0.; 
  
  
  //std vector
      ele_echarge.clear();
      //energy matches
      ele_he.clear(); 
      ele_hebc.clear(); 
      ele_eseedpout.clear(); 
      ele_ep.clear(); 
      ele_eseedp.clear(); 
      ele_eelepout.clear();       
      //delta 
      ele_deltaetaseed.clear(); 
      ele_deltaetaele.clear(); 
      ele_deltaphiseed.clear(); 
      ele_deltaphiele.clear(); 
      ele_deltaetain.clear(); 
      ele_deltaphiin.clear();
      //cluster shapes
      ele_sigmaietaieta.clear(); 
      ele_sigmaetaeta.clear(); 
      ele_sigmaiphiiphi.clear();; 
      ele_e15.clear(); 
      ele_e25max.clear(); 
      ele_e55.clear(); 
      ele_r9.clear();
      //momentum
      ele_pin_mode.clear(); 
      ele_pout_mode.clear(); 
      ele_pTin_mode.clear(); 
      ele_pTout_mode.clear(); 
      //brem
      ele_fbrem.clear();
      ele_SCfbrem.clear();
      ele_pfSCfbrem.clear(); 
      ele_trackfbrem.clear();
      ele_nbrem.clear();
      //position
      ele_isbarrel.clear(); 
      ele_isendcap.clear(); 
      ele_isEBetaGap.clear(); 
      ele_isEBphiGap.clear(); 
      ele_isEEdeeGap.clear(); 
      ele_isEEringGap.clear();
      ele_isecalDriven.clear(); 
      ele_istrackerDriven.clear();
      ele_eClass.clear();
      //distance vertex
      //ele_dxyB.clear(); 
      ele_dxy.clear(); 
      //ele_dzB.clear(); 
      ele_dz.clear(); 
      //ele_dszB.clear(); 
      ele_dsz.clear();              
      //conversion
      ele_valid_hits.clear(); 
      ele_lost_hits.clear(); 
      ele_chi2_hits.clear(); 
      ele_conv_dcot.clear();
      ele_conv_dist.clear();
      ele_conv_radius.clear();
      ele_expected_inner_hits.clear();
      //mva 
       ele_mva.clear();
      //isolation
      ele_pfChargedHadIso.clear(); 
      ele_pfNeutralHadIso.clear(); 
      ele_pfPhotonIso.clear(); 
      //supercluster energies
      ele_sclE.clear(); 
      ele_sclEt.clear(); 
      ele_sclEta.clear(); 
      ele_sclPhi.clear(); 
      ele_sclRawE.clear();
        ele_sclNclus.clear();
      //supercluster variables
       ele_sclphiwidth.clear(); 
       ele_scletawidth.clear();
      //errors
      ele_ecalErr.clear(); 
      ele_trackErr.clear(); 
      ele_combErr.clear(); 
      //kf       
       ele_kfchi2.clear();
       ele_kfhits.clear();
      //pf
      ele_pf_id.clear();
      ele_pf_eta.clear();
      ele_pf_phi.clear(); 
      ele_pf_pt.clear(); 
      ele_pf_dz.clear(); 
      ele_pf_dxy.clear();
      ele_pf_vx.clear(); 
      ele_pf_vy.clear(); 
      ele_pf_vz.clear();
      ele_pf_mva_nog.clear(); 
      ele_pf_mva_epi.clear();
      ele_pf_number.clear();
      
      ele_seed_subDet2.clear() ;
      ele_seed_dRz2.clear() ;
      ele_seed_dPhi2.clear() ;
      ele_seed_dRz2Pos.clear() ;
      ele_seed_dPhi2Pos.clear() ;
      ele_seed_subDet1.clear() ; 
      ele_seed_dRz1.clear();
      ele_seed_dPhi1.clear();
      ele_seed_dRz1Pos.clear();
      ele_seed_dPhi1Pos.clear();
  
      pho_pt.clear();
      pho_eta.clear();
      pho_phi.clear();
      pho_isEB.clear();
      pho_isEE.clear();
      pho_sigmaietaieta.clear();
      pho_sigmaetaeta.clear();
      pho_EcalRechitsSumEt.clear();  
      pho_HcalTowerSumEt.clear();
      pho_e15.clear();
      pho_e25.clear();
      pho_e33.clear();
      pho_e55.clear();
      pho_r9.clear();
      pho_hoe.clear();
      pho_sclEtaWidth.clear();
      pho_sclPhiWidth.clear();
      pho_sclRawE.clear();
      pho_pf_id.clear();
      pho_pf_eta.clear();
      pho_pf_phi.clear(); 
      pho_pf_pt.clear(); 
      pho_pf_dxy.clear();
      pho_pf_dz.clear();
      pho_pf_vx.clear();
      pho_pf_vy.clear();
      pho_pf_vz.clear();     
      pho_pf_mva_nog.clear(); 
      pho_pf_mva_epi.clear();
      pho_pf_number.clear();
      pho_isEBEEGap.clear(); 
      pho_isEBGap.clear(); 
      pho_isEBEtaGap.clear(); 
      pho_isEBPhiGap.clear();
      pho_isEEGap.clear(); 
      pho_isEEDeeGap.clear(); 
      pho_isEERingGap.clear();
      pho_isStandardPhoton.clear(); 
      pho_hasPixelSeed.clear();
      pho_maxEnergyXtal.clear();

}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Ntuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Ntuplizer);
