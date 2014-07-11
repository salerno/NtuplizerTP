import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntuplizer")

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")


process.GlobalTag.globaltag = 'START62_V1::All'
#START53_V15A::All'

process.load("Configuration.EventContent.EventContent_cff")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

# ----------------------------------------------------------------------
# Input File
# ----------------------------------------------------------------------
process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
    'file:000D3902-14E6-E311-8535-002354EF3BD0.root'
     #root://xrootd.unl.edu//store/relval/CMSSW_6_2_0_SLHC13_patch1/RelValZEE_14TeV/GEN-SIM-RECO/DES23_62_V1_UPG2023Muon-v1/00000/2AF00CAB-5DEA-E311-9793-0025905A6060.root'
                )
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
                                   VerticesTag = cms.InputTag('offlinePrimaryVertices'),
                                   TracksTag = cms.InputTag('generalTracks'),
                                   isMC = cms.bool(True)
                                   )

# ----------------------------------------------------------------------
# Path to execute
# ----------------------------------------------------------------------
process.p = cms.Path(process.ntuplizer)
