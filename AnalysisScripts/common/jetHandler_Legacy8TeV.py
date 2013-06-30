import FWCore.ParameterSet.Config as cms

## PU Jet ID
from CMGTools.External.pujetidproducer_cfi import simple, full, cutbased

### ## JEC
jec = cms.PSet(
    data = cms.vstring("aux/jec2012/FT_53_V21_AN4::All_L1FastJet_AK5PF.txt",
                       "aux/jec2012/FT_53_V21_AN4::All_L2Relative_AK5PF.txt",
                       "aux/jec2012/FT_53_V21_AN4::All_L3Absolute_AK5PF.txt",
                       "aux/jec2012/FT_53_V21_AN4::All_L2L3Residual_AK5PF.txt",
                       ),
    mc   = cms.vstring("aux/jec2012/START53_V23::All_L1FastJet_AK5PF.txt",
                       "aux/jec2012/START53_V23::All_L2Relative_AK5PF.txt",
                       "aux/jec2012/START53_V23::All_L3Absolute_AK5PF.txt",
                       ),
    unc = cms.string("aux/jec2012/START53_V23::All_Uncertainty_AK5PF.txt")
    )

jer = cms.PSet(
    etaBins  = cms.vdouble(0,0.5,1.1,1.7,2.3),
    resSf    = cms.vdouble(1.052,1.057,1.096,1.134,1.288),
    resSfErr = cms.vdouble(0.063,0.057,0.065,0.093,0.200)
    )
