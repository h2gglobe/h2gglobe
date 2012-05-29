import FWCore.ParameterSet.Config as cms

## PU Jet ID
from CMGTools.External.pujetidproducer_cfi import simple, full, cutbased

### ## JEC
jec = cms.PSet(
    data = cms.vstring("aux/jec2012/Summer12_V1_DATA_L1FastJet_AK5PF.txt",
                       "aux/jec2012/Summer12_V1_DATA_L2Relative_AK5PF.txt",
                       "aux/jec2012/Summer12_V1_DATA_L3Absolute_AK5PF.txt",
                       "aux/jec2012/Summer12_V1_DATA_L2L3Residual_AK5PF.txt",
                       ),
    mc   = cms.vstring("aux/jec2012/Summer12_V1_MC_L1FastJet_AK5PF.txt",
                       "aux/jec2012/Summer12_V1_MC_L2Relative_AK5PF.txt",
                       "aux/jec2012/Summer12_V1_MC_L3Absolute_AK5PF.txt",
                       )
    )
