import FWCore.ParameterSet.Config as cms

DQMExample_Step1 = cms.EDAnalyzer('DemoAnalyzer', 
     processFastSim = cms.untracked.bool(False),
     SimHitTagLabels = cms.VInputTag(
        #cms.InputTag("famosSimHits", "TrackerHits")
        cms.InputTag("g4SimHits", "TrackerHitsPixelBarrelLowTof"),
        cms.InputTag("g4SimHits", "TrackerHitsPixelBarrelHighTof"),
        cms.InputTag("g4SimHits", "TrackerHitsPixelEndcapLowTof"), 
        cms.InputTag("g4SimHits", "TrackerHitsPixelEndcapHighTof"),
        cms.InputTag("g4SimHits", "TrackerHitsTECLowTof"), 
        cms.InputTag("g4SimHits", "TrackerHitsTECHighTof"),
        cms.InputTag("g4SimHits", "TrackerHitsTIBLowTof"), 
        cms.InputTag("g4SimHits", "TrackerHitsTIBHighTof"),
        cms.InputTag("g4SimHits", "TrackerHitsTIDLowTof"), 
        cms.InputTag("g4SimHits", "TrackerHitsTIDHighTof"),
        cms.InputTag("g4SimHits", "TrackerHitsTOBLowTof"), 
        cms.InputTag("g4SimHits", "TrackerHitsTOBHighTof")
        ),
     bins_p = cms.untracked.vstring('2.1','1.2', '3.4', '.8', '5.2', '.4'),
    )



