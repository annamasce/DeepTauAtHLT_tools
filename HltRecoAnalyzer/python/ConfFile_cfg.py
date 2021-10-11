import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            'file:/afs/cern.ch/work/a/amascell/CMSSW_11_2_0/src/myFile_inMINIAODSIM.root'
                )
                            )

process.demo = cms.EDAnalyzer('HltRecoAnalyzer',
    hltresults = cms.InputTag("TriggerResults", "", "reRECO"),
    trig_objects = cms.InputTag("slimmedPatTrigger", "", "reRECO"),
    taus_offline = cms.InputTag("slimmedTaus"),
    gen_particles = cms.InputTag("genParticles"),
    hlt_path = cms.string("HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4"),
    filter_to_pass = cms.string("hltHpsDoublePFTau30TrackPt1DeepTau30IsolationDz02Reg"),
    run_parameters = cms.PSet(
        HistogramFile = cms.untracked.string("HltRecoInfo.root")
    )
)

process.p = cms.Path(process.demo)