import FWCore.ParameterSet.Config as cms

process = cms.Process("Validation")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.EDMtoMEAtRunEnd_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

f = open('muons_FULLSIM_GENSIM_02_5.txt', 'r')

myfilelist = cms.untracked.vstring()
myfilelist.extend( [line.strip() for line in f.read().splitlines()] )
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames =
        #myfilelist
        cms.untracked.vstring(
        #FullSim Samples (locally)
        #'file:/afs/cern.ch/user/g/gkrintir/github/GenSim/CMSSW_7_3_0/src/mygun_ppions_FULLSIM_GENSIM.root'
        
        #FullSim Samples (on the grid)
        
        #FastSim Samples (locally)
        #'file:/afs/cern.ch/user/g/gkrintir/github/GenSim/CMSSW_7_3_0/src/mygun_muons_FASTSIM_SameFULLSIMConditions.root'

        'file:/afs/cern.ch/user/g/gkrintir/github/GenSim/CMSSW_7_3_0/src/Generation_Output/FastSim/GENSIM/mygun_ppions_FASTSIM_SameFULLSIMConditions_noflatpT.root',
        #'file:/afs/cern.ch/user/g/gkrintir/github/GenSim/CMSSW_7_3_0/src/mygun_protons_FASTSIM_SameFULLSIMConditions_noflatpT.root,
        #'file:/afs/cern.ch/user/g/gkrintir/github/GenSim/CMSSW_7_3_0/src/mygun_kaons_FASTSIM_SameFULLSIMConditions_noflatpT.root'
        )
)

# For histograms
process.load("DQMServices.Core.DQM_cfg")
process.load("DQMServices.Components.DQMEnvironment_cfi")

process.validation = cms.EDAnalyzer('DemoAnalyzer', 
     processFastSim = cms.untracked.bool(True),                               
)

process.GlobalTag.globaltag = 'PHYS14_50_V1'

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('nalysis_mygyn_muons_FASTSIM_SameFULLSIMConditions_noflatpT.root')
                                   )

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
                                     fileName = cms.untracked.string("OUT_step1.root"))

#process.p = cms.Path(process.validation)

# Path and EndPath definitions
process.dqmoffline_step = cms.Path(process.validation)
'''
process.dqmsave_step = cms.Path(process.DQMSaver)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)


# Schedule definition
process.schedule = cms.Schedule(
    process.dqmoffline_step,
    process.DQMoutput_step,
    process.dqmsave_step
    )
'''


process.dqmsave_step = cms.Path(process.DQMSaver)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)


# Schedule definition
process.schedule = cms.Schedule(
    process.dqmoffline_step,
    process.DQMoutput_step,
    process.dqmsave_step
    )
