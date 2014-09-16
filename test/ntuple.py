import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntuplizer")

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryExtended2023MuonReco_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.GlobalTag.globaltag = 'START62_V1::All'
#START53_V15A::All'

process.load("Configuration.EventContent.EventContent_cff")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))

# ----------------------------------------------------------------------
# Input File
# ----------------------------------------------------------------------
process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
	
fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/mc/GEM2019Upg14DR/GluGluToHToGG_M-125_14TeV-powheg-pythia6/AODSIM/final_phase1_PU50bx25_DES19_62_V8-v1/30000/007E9600-9E1B-E411-8E36-0025905B85EE.root')	
#fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/mc/GEM2019Upg14DR/GJet_Pt-15to3000_Tune4C_14TeV_pythia8/AODSIM/final_phase1_age1k_PU140bx25_PH1_1K_FB_V2-v1/00000/00493C8C-8E1F-E411-B826-0025905A60CE.root')	                             
    #fileNames = cms.untracked.vstring(
    #'file:drellyann.root'
    #'file:10AF42EE-1AE2-E311-8F1E-00261894382A.root'
    #'file:RelVal/relval.root'
    #'file:GJets.root'
     #root://xrootd.unl.edu//store/relval/CMSSW_6_2_0_SLHC13_patch1/RelValZEE_14TeV/GEN-SIM-RECO/DES23_62_V1_UPG2023Muon-v1/00000/2AF00CAB-5DEA-E311-9793-0025905A6060.root'
     #          )
                            )

# ----------------------------------------------------------------------
# Output root file (monitoring histograms)
# ----------------------------------------------------------------------
process.TFileService=cms.Service('TFileService',
                                fileName=cms.string('TP_ntuple.root')
                                )


# ----------------------------------------------------------------------
# Ntuplizer
# ----------------------------------------------------------------------
process.ntuplizer = cms.EDAnalyzer('Ntuplizer',
                                   EleTag      = cms.InputTag('gsfElectrons'),
                                   PhoTag      = cms.InputTag('photons'),
				   BarrelSCTag = cms.InputTag('particleFlowSuperClusterECAL:particleFlowSuperClusterECALBarrel'),
				   EndcapsSCTag = cms.InputTag('particleFlowSuperClusterECAL:particleFlowSuperClusterECALEndcapWithPreshower'),
                                   VerticesTag = cms.InputTag('offlinePrimaryVertices'),
                                   TracksTag = cms.InputTag('generalTracks'),
                                   isMC = cms.bool(True)
                                   )

# ----------------------------------------------------------------------
# Path to execute
# ----------------------------------------------------------------------
process.p = cms.Path(process.ntuplizer)
