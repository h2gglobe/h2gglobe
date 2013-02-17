templates = (
    'pho{pho}_sigmaEOverE_{ebee2}{syst}_cat0_DYJetsToLL',
    'pho{pho}_sigmaEOverE_{R9}R9{syst}_cat{cat8}_DYJetsToLL',
    'pho{pho}_phoidMva_{ebee2}{syst}_cat0_DYJetsToLL',
    'pho{pho}_phoidMva_{R9}R9{syst}_cat{cat8}_DYJetsToLL',
    'bdtout{ebee}{syst}_cat0_DYJetsToLL',
    'bdtout{syst}_cat{cat4}_DYJetsToLL',
    'bdtout_{R9t}R9{syst}_cat{cat4}_DYJetsToLL',
    )

# 'syst' is a special variation and should correspond to the things to be morphed: nominal, plus, minus
variations = {
    'pho' : (1,2),
    'R9' : ('low','mid','high'),
    'cat8' : tuple(range(8)),
    'cat4' : tuple(range(4)),
    'ebee' : ('EB','EE','EBEE'),
    'ebee2' : ('EB','EE'),
    'R9t' : ('low','mixed','high'),
    'syst' : ('','_up','_down'),
    }
