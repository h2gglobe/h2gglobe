import FWCore.ParameterSet.Config as cms

## PU Jet ID
from CMGTools.External.pujetidproducer_cfi import simple, full, cutbased

### ## JEC
jec = cms.PSet(
    data = cms.vstring("aux/jec2012/Summer12_V1_DATA::All_L1FastJet_AK5PF_new.txt",
                       "aux/jec2012/Summer12_V1_DATA::All_L2Relative_AK5PF_new.txt",
                       "aux/jec2012/Summer12_V1_DATA::All_L3Absolute_AK5PF_new.txt",
                       "aux/jec2012/Summer12_V1_DATA::All_L2L3Residual_AK5PF.txt",
                       ),
    mc   = cms.vstring("aux/jec2012/Summer12_V1_MC::All_L1FastJet_AK5PF_new.txt",
                       "aux/jec2012/Summer12_V1_MC::All_L2Relative_AK5PF_new.txt",
                       "aux/jec2012/Summer12_V1_MC::All_L3Absolute_AK5PF_new.txt",
                       )
    )
