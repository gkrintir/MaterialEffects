import FWCore.ParameterSet.Config as cms

DQMExample_Step1 = cms.EDAnalyzer('DemoAnalyzer', 
     processFastSim = cms.untracked.bool(False),
     SimHitTags = cms.VInputTag(
        cms.InputTag("famosSimHits", "TrackerHits")
        #cms.InputTag("g4SimHits", "TrackerHitsPixelBarrelLowTof"),
        #cms.InputTag("g4SimHits", "TrackerHitsPixelBarrelHighTof"),
        #cms.InputTag("g4SimHits", "TrackerHitsPixelEndcapLowTof"), 
        #cms.InputTag("g4SimHits", "TrackerHitsPixelEndcapHighTof"),
        #cms.InputTag("g4SimHits", "TrackerHitsTECLowTof"), 
        #cms.InputTag("g4SimHits", "TrackerHitsTECHighTof"),
        #cms.InputTag("g4SimHits", "TrackerHitsTIBLowTof"), 
        #cms.InputTag("g4SimHits", "TrackerHitsTIBHighTof"),
        #cms.InputTag("g4SimHits", "TrackerHitsTIDLowTof"), 
        #cms.InputTag("g4SimHits", "TrackerHitsTIDHighTof"),
        #cms.InputTag("g4SimHits", "TrackerHitsTOBLowTof"), 
        #cms.InputTag("g4SimHits", "TrackerHitsTOBHighTof")
     ),
     #particles = cms.untracked.string('pions/'), #token(e.g. :) could be everything other a number, a letter and an underscore
     particleTypes = cms.untracked.vint32(211),
     #bins_p =  cms.untracked.vdouble(2.1, 1.2, 3.4, .8, 5.2, .4),
     bins_p =  cms.untracked.vdouble( 1.2, .4, 1.5, 1.9),
     TestParticleFilter = cms.PSet(
        # Particles with |eta| > etaMax (momentum direction at primary vertex)
        # are not simulated 
        etaMin = cms.double(-2.5),
        # Particles with |eta| > etaMax (momentum direction at primary vertex)
        # are not simulated
        etaMax = cms.double(2.5),
        # Charged particles with pT < pTMin (GeV/c) are not simulated
        pTMin = cms.double(0.05),
        # Charged particles with pT < pTMin (GeV/c) are not simulated
        pTMax = cms.double(100.0),
        # Charged particles with pT < pTMin (GeV/c) are not simulated
        pMin = cms.double(0.0),
        # Charged particles with pT < pTMin (GeV/c) are not simulated
        pMax = cms.double(0.0),
        # Particles with energy smaller than EMin (GeV) are not simulated
        EMin = cms.double(0.0),
        # Particles with energy smaller than EMin (GeV) are not simulated
        EMax = cms.double(0.0),
        # Protons with energy in excess of this value (GeV) will kept no matter what
        pdgIdsToFilter = cms.vint32(2212)
        ),
     formating1D = cms.VPSet (
        cms.PSet(
         title = cms.string('mine'),
         name = cms.string('default'),
         labelx = cms.untracked.string('dEdX'),
         labely = cms.untracked.string(''),
         rangex = cms.untracked.vdouble(200, 0., 1e-2)
         ),
        ),
     formating2D = cms.VPSet (
        cms.PSet(
         title = cms.string('default'),
         name = cms.string('default'),
         labelx = cms.untracked.string('p (#it{GeV}/c)'),
         labely = cms.untracked.string('dEdX'),
         rangex = cms.untracked.vdouble(60.,120.),
         rangey = cms.untracked.vdouble(60.,120.),
         )
        ) 
 )                          

    



