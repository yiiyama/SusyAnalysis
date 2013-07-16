import FWCore.ParameterSet.Config as cms

softLoosePhotonFilter = cms.EDFilter("SoftLoosePhotonFilter",
    photonCollectionTag = cms.untracked.InputTag("photons"),
    electronCollectionTag = cms.untracked.InputTag("gsfElectrons"),
    conversionCollectionTag = cms.untracked.InputTag("allConversions"),
    beamSpotTag = cms.untracked.InputTag("offlineBeamSpot"),
    ptThreshold = cms.untracked.double(15.),
    cuts = cms.untracked.vint32()
)
