#!/usr/bin/env python

import sys, os, re

#----------------------------------------------------------------------

def floatEquals(x,y):
    return abs(x-y) < 1e-4

#----------------------------------------------------------------------
# main
#----------------------------------------------------------------------
from optparse import OptionParser

parser = OptionParser("""

%prog [options] input_file1.root [ input_file2.root ... ]
%prog [options] input_file1.log [ input_file2.log ... ]

    takes the observed and expected limit from the
    per-mass point ROOT files (output of CMS' combine program)
    and produces a CSV file (by default written to stdout)
"""
)

parser.add_option(
    "--output",
    dest="outputFname",
    help="instead of writing to stdout, write to the given file. ",
    default=None,
    type = str,
    )

parser.add_option(
    "--ignore-missing-numbers",
    dest="ignoreMissingNumbers",
    help="do not return with a non-zero exit value if some of the numbers were missing (e.g. the median expected limit)",
    default=False,
    action = "store_true",
    )

parser.add_option(
    "--clb",
    help="read CLb values from the log files (i.e. log files must be specified, not root files).",
    default=False,
    action = "store_true",
    )
(options, ARGV) = parser.parse_args()

#----------------------------------------

if options.outputFname != None:
    if os.path.exists(options.outputFname):
        print >> sys.stderr,"output file '%s' already exists, refusing to overwrite it" % options.outputFname
        sys.exit(1)

    outputFile = open(options.outputFname,"w")
else:
    outputFile = sys.stdout

# print the header
if options.clb:
    header = ['mass','CLb']
else:
    header = ['mass','xsect_lim_obs','xsect_lim_expected',
    
                    'xsect_lim_expected_plus_1_sigma',
                    'xsect_lim_expected_minus_1_sigma',
    
                    'xsect_lim_expected_plus_2_sigma',
                    'xsect_lim_expected_minus_2_sigma',
                    ]

print >> outputFile, ",".join(header)

import ROOT

errFlag = False

for fname in ARGV:
    
    if options.clb:
        clb = None

        row = [ -1 ] * len(header)
        
        # extract the mass from the file name
        fnameLastPart = fname[fname.rfind('/')+1:]
        
        mo = re.search("([0-9]+\.[0-9])", fname)
        assert(mo)
        
        mass = float(mo.group(1))
        
        # read clb from the log file
        lines = open(fname).readlines()
        
        # example line:
        # At r = 1.000000:   q_mu = 0.71362  q_A  = 1.95085  CLsb = 0.19912  CLb  = 0.70951  CLs  = 0.28065
    
        # take the LAST line, i.e. start looking at the lines from the end
        for line in lines[::-1]:
            mo = re.match("At r = \S+:\s+q_mu = \S+\s+q_A\s+= \S+\s+CLsb = \S+\s+CLb\s+= (\S+)\s+CLs\s+= \S+$",line)
            
            if mo:
                clb = mo.group(1)
            
            
        # insist that we found such a line
        assert(clb != None)
        row[header.index('CLb')] = clb
    
    else:
        # read the exclusion limits 
        # mo = re.match("higgsCombineTest.Asymptotic.mH(.+).root",fname)
        # assert(mo)
        # mass = float(mo.group(1))
    
        fin = ROOT.TFile(fname)
    
        tree = fin.Get("limit")
        tree.Draw("limit:mh:quantileExpected","","goff")
    
        numSelectedRows = tree.GetSelectedRows()
    
        row = [ -1 ] * len(header)
        mass = tree.GetV2()[0]
    
        # assert(numSelectedRows == 6)
        if numSelectedRows != 6:
            print >> sys.stderr,"warning: not all limit values (expected/observed) found for mass %.1f in file %s" % (mass, fname)
            if not options.ignoreMissingNumbers:
                errFlag = True
    
    
        for i in range(numSelectedRows):
    
            quantile = tree.GetV3()[i]
            limit = tree.GetV1()[i]
    
    
            if floatEquals(quantile, -1):
                row[header.index('xsect_lim_obs')] = limit
            elif floatEquals(quantile, 0.025):
                row[header.index('xsect_lim_expected_minus_2_sigma')] = limit
            elif floatEquals(quantile, 0.16):
                row[header.index('xsect_lim_expected_minus_1_sigma')] = limit
            elif floatEquals(quantile, 0.5):
                # median expected
                row[header.index('xsect_lim_expected')] = limit
    
            elif floatEquals(quantile, 0.84):
                row[header.index('xsect_lim_expected_plus_1_sigma')] = limit
            elif floatEquals(quantile, 0.975):
                row[header.index('xsect_lim_expected_plus_2_sigma')] = limit
            else:
                raise Exception("don't know quantile" + str(quantile))
    
        # OLD:
        # indices were found empirically
        # by comparing ntuple content and program output
        # observedRatioLimit = tree.GetV1()[5]
        # expectedRatioLimit = tree.GetV1()[2]

    row[header.index('mass')] = mass
    

    # write the values to the output file
    print >> outputFile, ",".join([str(x) for x in row ])
    
    
if options.outputFname != None:
    outputFile.close()


if errFlag:
    # an error occurred
    sys.exit(1)
