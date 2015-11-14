import FWCore.ParameterSet.Config as cms

process = cms.Process('RECODQM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EDMtoMEAtRunEnd_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

# load DQM
process.load("DQMServices.Core.DQM_cfg")
process.load("DQMServices.Components.DQMEnvironment_cfi")
#process.load("DQMServices.Components.EDMtoMEConverter_cff")
# my analyzer
process.load("Analyzer.MaterialEffectsAnalyzer.DQMExample_Step1_cfi")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(500)
)

f = open('ppions_FULLSIM_GENSIM_05_10.txt', 'r')

myfilelist = cms.untracked.vstring()
myfilelist.extend( [line.strip() for line in f.read().splitlines()] )

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    #fileNames = myfilelist
    fileNames = cms.untracked.vstring(
        #reco from relVals
        #'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/datatest/forTutorial/step2_RAW2DIGI_RECO_fromRelValTTbarLepton.root'
        'file:/afs/cern.ch/user/g/gkrintir/github/GenSim/CMSSW_7_3_0/src/Generation_Output/FastSim/GENSIM/mygun_ppions_FASTSIM_SameFULLSIMConditions_noflatpT.root'
        )
)

process.GlobalTag.globaltag = 'PHYS14_50_V1'

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('nalysis_mygyn_muons_FASTSIM_SameFULLSIMConditions_noflatpT.root')
                                   )

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
                                     fileName = cms.untracked.string("OUT_step1.root"))
'''
# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')  #for MC
'''

#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup_GRun', '')

process.dqmSaver.workflow = "/CMSSW_3_1_1/RelVal/TrigVal_Hgg"
#process.dqmSaver.convention = 'Offline'
# Path and EndPath definitions
process.dqmoffline_step = cms.Path(process.DQMExample_Step1)
process.dqmsave_step = cms.Path(process.dqmSaver)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)


# Schedule definition
process.schedule = cms.Schedule(
    process.dqmoffline_step,
    process.DQMoutput_step,
    process.dqmsave_step
    )

# Keep the logging output to a nice level
process.MessageLogger.destinations = ['newE_.txt']
