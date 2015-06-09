import FWCore.ParameterSet.Config as cms

process = cms.Process("Validation")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

f = open('kaons_FULLSIM_GENSIM_01_5.txt', 'r')

myfilelist = cms.untracked.vstring()
myfilelist.extend( [line.strip() for line in f.read().splitlines()] )
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames =
        myfilelist
        #cms.untracked.vstring(
        #FullSim Samples (locally)
        #'file:/afs/cern.ch/user/g/gkrintir/github/GenSim/CMSSW_7_3_0/src/mygun_ppions_FULLSIM_GENSIM.root'
        
        #FullSim Samples (on the grid)
        
        #FastSim Samples (locally)
        #'file:/afs/cern.ch/user/g/gkrintir/github/GenSim/CMSSW_7_3_0/src/mygun_ppions_FASTSIM_SameFULLSIMConditions.root'

        #'file:/afs/cern.ch/user/g/gkrintir/github/GenSim/CMSSW_7_3_0/src/Generation_Output/FastSim/GENSIM/mygun_kaons_FASTSIM_SameFULLSIMConditions_noflatpT.root',
        #'file:/afs/cern.ch/user/g/gkrintir/github/GenSim/CMSSW_7_3_0/src/mygun_protons_FASTSIM_SameFULLSIMConditions_noflatpT.root,
        #'file:/afs/cern.ch/user/g/gkrintir/github/GenSim/CMSSW_7_3_0/src/mygun_kaons_FASTSIM_SameFULLSIMConditions_noflatpT.root'
        #)
)

process.validation = cms.EDAnalyzer('DemoAnalyzer', 
     processFastSim = cms.untracked.bool(False),                               
)

process.GlobalTag.globaltag = 'PHYS14_50_V1'

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('Analysis_mygyn_kaons_FASTSIM_SameFULLSIMConditions_noflatpT.root')
                                   )


process.p = cms.Path(process.validation)
