#!/usr/bin/env python

import time
import os
import csv

import PlotUtils

import sys
#----------------------------------------------------------------------
# parameters
#----------------------------------------------------------------------

colorsToAvoid = [ 5, # yellow
                  ] 


#----------------------------------------------------------------------

def makePlot(csvFnames, relative, includeExpected = True, fermiophobic = None, ymax = None, inputIsAbs = None, drawXsectBR = False,
             minMass = None,
             maxMass = None,
             plotLog = False
             ):
    """ @param relative if true, the exclusion of the ratios
        (relative to the inputs given) are plotted. If False,
        these ratios are multiplied by the SM cross sections
        and branching ratios into gamma gamma 

        @param inputIsAbs is a list (set) of file names which
        should be treated as if they had ABSOLUTE limits on
        cross sections rather than relative limits.

        @param minMass and masMass can be used to restrict the
        plotting range
    """

    #--------------------
    # read the files            
    #--------------------
    data = []
    color = 0

    for fname in csvFnames:

        while True:
            color += 1
            if color not in colorsToAvoid:
                break

        # define a name: if there is a comma in the file name assume that the
        # part before the comma is the actual file name and the part after it
        # is the label we should use
        #
        # if there is no comma, just use the basename (without .csv) as label
        pos = fname.find(',')
        if pos == -1:
            # not found
            name = os.path.basename(fname)
            if name.lower().endswith(".csv"):
                name = name[:-4]
        else:
            name = fname[pos+1:]
            fname = fname[:pos]

        masses, observedValues, expectedValues, \
           expected_minus_2_sigma_values, \
           expected_minus_1_sigma_values, \
           expected_plus_1_sigma_values, \
           expected_plus_2_sigma_values =  PlotUtils.readCSV(open(fname), includeExpected)

        #--------------------
        # filter on masses
        #--------------------
        indices = range(len(masses))

        if minMass != None:
            indices = [ i for i in indices if masses[i] >= minMass ]

        if maxMass != None:
            indices = [ i for i in indices if masses[i] <= maxMass ]

        masses = [ masses[i] for i in indices ]
        observedValues = [ observedValues[i] for i in indices ]
        expectedValues = [ expectedValues[i] for i in indices ]

        #--------------------
        
        tmp = { "masses": masses,
                "observedValues": observedValues,
                "expectedValues": expectedValues,

                # for labels
                "name": name,

                # assign the color here
                "color": color,
                }


            
        data.append(tmp)    


    #--------------------
            
    # just to make sure we're not picking up something in the code afterwards
    del masses
    del observedValues
    del expectedValues
    #--------------------
            
    if not relative:

        # if we're plotting the absolute cross sections, we
        # need to know whether this is Standard Model or Fermiophobic
        assert(fermiophobic != None)

        if fermiophobic:
            typeName = "FP"
        else:
            typeName = "SM"

        # convert to absolute cross sections
        for line, fname in zip(data, csvFnames):
            if fname in inputIsAbs:
                # already absolute
                continue

            line['observedValues'] = PlotUtils.multiplyArrayByXsectAndBR(line['masses'], line['observedValues'], fermiophobic)
            line['expectedValues'] = PlotUtils.multiplyArrayByXsectAndBR(line['masses'], line['expectedValues'], fermiophobic)

    else:
        # we're asked to plot relative results, convert to relative for those
        # inputs which are absolute

        for line, fname in zip(data, csvFnames):
            if not fname in inputIsAbs:
                # relative input, no need to convert
                continue

            line['observedValues'] = PlotUtils.divideArrayByXsectAndBR(line['masses'], line['observedValues'], fermiophobic)
            line['expectedValues'] = PlotUtils.divideArrayByXsectAndBR(line['masses'], line['expectedValues'], fermiophobic)

    #----------------------------------------
    # legend
    #----------------------------------------
    legend = ROOT.TLegend(options.legendXleft, options.legendYbottom,
                          options.legendXright,options.legendYtop); gcSaver.append(legend)
    
    legend.SetShadowColor(0);
    legend.SetFillColor(0);
    legend.SetBorderSize(0);

    #----------------------------------------
    # produce the 'observed' graphs
    #----------------------------------------

    allGraphs = []

    for line in data:
        gr = PlotUtils.makeGraphFromArrays(line['masses'], line['observedValues'])
        line['grObserved'] = gr
        gcSaver.append(gr)

        if options.observedLineWidth > 0:
            gr.SetLineWidth(options.observedLineWidth)
        else:
            # set default width for legend
            gr.SetLineWidth(4)

        gr.SetLineColor(line['color'])

        legend.AddEntry(gr,line['name'],"L")

        if options.observedLineWidth > 0:
            allGraphs.append(gr)

    #----------------------------------------
    # produce the 'expected' graphs
    #----------------------------------------

    if includeExpected:

        for line in data:
            
            grExpected = PlotUtils.makeGraphFromArrays(line['masses'], line['expectedValues'])
            gcSaver.append(grExpected)

            line['grExpected'] = grExpected

            grExpected.SetLineStyle(ROOT.kDashed)
            grExpected.SetLineWidth(4)
            grExpected.SetLineColor(line['color'])

            allGraphs.append(grExpected)

            # label = makeGraphLabelOnRight(grExpected, minMass, maxMass, "BG exp.")
            # label.SetTextSize(label.GetTextSize() * 0.7)    
            # label.Draw()
            # gcSaver.append(label)

    #myCanvas = ROOT.TCanvas("myCanvas","Title Goes Here")
    #myCanvas.SetLogy(plotLog)

    #----------------------------------------
    # produce the graph for the theoretical cross section
    #----------------------------------------
    if drawXsectBR:
        # add a graph for the theoretical cross section

        # take the 'or' of all masses given
        import operator
        allMasses = sorted(reduce(operator.add, [ line['masses'] for line in data ] ))

        # for the moment, limit this to integer masses (in GeV)
        # (the cross section interpolation seems not yet to be functional)
        allMasses = sorted(list(set([ int(round(x)) for x in allMasses ])))

        # print "allMasses=",allMasses

        theoXsectBr = [ PlotUtils.getXsectTimesBR(mass, fermiophobic) for mass in allMasses ]
        gr = PlotUtils.makeGraphFromArrays(allMasses, theoXsectBr)

        gr.SetLineWidth(4)
        gr.SetLineStyle(ROOT.kDotted)

        legend.AddEntry(gr,"theo. #sigma * BR","L")

        gcSaver.append(gr)
        allGraphs.append(gr)

    #----------------------------------------
    # determine the y scale
    #----------------------------------------

    if ymax == None:

        # determine this from the values, not from the graphs
        # (is there a way to do this from the graphs ?)
        ymax = max([value for line in data for value in line['observedValues'] ]) 

        if includeExpected:
            ymax = max(ymax, max([value for line in data for value in line['expectedValues'] ]))

        ymax *= 1.1

        # TODO: remove this if statement ?!
        if not relative:
            if fermiophobic:
                # fix the y scale by hand in order not to
                # stretch it too much because of large
                # scaling factors for the theoretical expectation
                ymax = 0.5

    #----------------------------------------
    # determine x scale (mass range)
    #----------------------------------------
    allMasses = [value for line in data for value in line['masses'] ]
    actualMinMass = min(allMasses)
    actualMaxMass = max(allMasses)

    del allMasses

    #----------------------------------------

    # create a dummy histogram to set the x range
    hdummy = ROOT.TH1F("hdummy","",1,actualMinMass,actualMaxMass)
    gcSaver.append(hdummy)

    hdummy.SetMaximum(ymax)
    hdummy.Draw()

    ROOT.gStyle.SetOptTitle(0)

    #----------------------------------------
    # draw the graphs
    #----------------------------------------

    for gr in allGraphs:
        gr.Draw("C,same")
        #gr.Draw("L,same")

        
    #----------------------------------------

    ROOT.gStyle.SetOptStat(0)
    hdummy.SetXTitle("m_{H} [GeV/c^{2}]")
    hdummy.GetYaxis().SetTitleOffset(1.2 * hdummy.GetYaxis().GetTitleOffset())

    if relative:
        hdummy.SetYTitle("#sigma/#sigma(theo)")
    else:
        hdummy.SetYTitle("#sigma(%s) * BR(#gamma#gamma) [pb]" % typeName)
        
    ROOT.gPad.SetGrid()

    if options.showTitle:
        label = ROOT.TLatex(0.5,0.85,"Excluded at 95% CL.")
        gcSaver.append(label)
        label.SetNDC(1)

        label.SetTextAlign(21)
        label.Draw()


    legend.Draw()

    ROOT.gPad.Modified()

    ROOT.gPad.Modified()


                    

#----------------------------------------------------------------------
# main
#----------------------------------------------------------------------
from optparse import OptionParser

parser = OptionParser("""

%prog [options] csv_file1 csv_file2 [ ... ]

    compares two or more cross section ratio limit
    output files on one plot
"""
)

parser.add_option(
    "--fermiophobic",
    dest="fermiophobicMode",
    help="exclude gluon fusion and ttbar associated production and rescale the other two processes such that " + \
      "their total cross section corresponds to the previous total",
    default=False,
    action="store_true",
    )

parser.add_option(
    "--save",
    dest="outputFnames",
    help="save the plot into a file (can be specified multiple times)",
    default=[],
    action="append",
    )

parser.add_option(
    "--plotLog",
    help="plots y in log scale",
    default=False,
    action="store_true",
    )

parser.add_option(
    "--ymax",
    type=float,
    help="manually sepcify the y scale",
    default=None,
    )

parser.add_option(
    "--relative",
    help="instead of plotting the absolute cross section exclusions, plot the relative (w.r.t. to the input signal)",
    default=False,
    action="store_true",
    )

parser.add_option(
    "--isabs",
    help="specify that a given file contains ABSOLUTE rather than RELATIVE (w.r.t. to the standard cross section) limits. Files can either be specified by name (which must be the same string as given in the list of non-option arguments) or by position (starting from 1)",
    default=[],
    action="append",
    )

parser.add_option(
    "--theo-xsect",
    help="add a line for the theoretical cross section times branching ratio",
    default=False,
    action="store_true",
    dest = "theoXsect",
    )

parser.add_option(
    "--min-mass",
    help="do not include values below the given mass",
    default=None,
    type=float,
    dest = "minMass",
    )

parser.add_option(
    "--max-mass",
    help="do not include values above the given mass",
    default=None,
    type=float,
    dest = "maxMass",
    )

parser.add_option(
    "--legend-xleft",
    help="NDC position of left side of legend",
    default=0.7,
    type=float,
    dest = "legendXleft",
    )

parser.add_option(
    "--legend-xright",
    help="NDC position of right side of legend",
    default=0.9,
    type=float,
    dest = "legendXright",
    )


parser.add_option(
    "--legend-ybottom",
    help="NDC position of bottom side of legend",
    default=0.7,
    type=float,
    dest = "legendYbottom",
    )

parser.add_option(
    "--legend-ytop",
    help="NDC position of top side of legend",
    default=0.9,
    type=float,
    dest = "legendYtop",
    )

parser.add_option(
    "--no-title",
    help="disable plotting of the title",

    # note the inverted logic here
    default = True,
    dest = "showTitle",
    action = "store_false", 
    )


parser.add_option(
    "--observed-line-width",
    help="line width for observed graphs (set to zero to NOT ot show them)",
    type = int,
    default = 4,
    dest = "observedLineWidth",
    )

(options, ARGV) = parser.parse_args()

#----------------------------------------
# process command line arguments
#----------------------------------------

if len(ARGV) < 1:
    PlotUtils.optError("expected at least one non-option arguments.")

# check whether any files were specified to contain ABSOLUTE cross
# sections
isabs = options.isabs
options.isabs = set()

for value in isabs:
    # check whether this is a file name
    if value in ARGV:
        # arriving here, it's a file name we know
        options.isabs.add(value)
        continue

    # try by position
    try:
        index = int(value)-1
        options.isabs.add(ARGV[index])
        continue
    except ValueError:
        # string did not represent a valid integer
        pass
    except IndexError:
        # out of index in array access
        pass

    print >> sys.stderr,"--isabs argument '%s' is neither a filename given nor a valid position" % value
    sys.exit(1)


#----------------------------------------------------------------------

import ROOT
ROOT.gROOT.SetStyle("Plain")
gcSaver = []

# makePlot(open(ARGV[0]), relative = True)
# ROOT.gPad.SaveAs("xsect-ratio-exclusion.eps")
# ROOT.gPad.SaveAs("xsect-ratio-exclusion.png")

inputFnames = ARGV[:]

makePlot(inputFnames, relative = options.relative, fermiophobic = options.fermiophobicMode, ymax = options.ymax, inputIsAbs = options.isabs,
         drawXsectBR = options.theoXsect,
         minMass = options.minMass,
         maxMass = options.maxMass,
         plotLog = options.plotLog,
         )

if options.outputFnames:
    for fname in options.outputFnames:
        print fname
        ROOT.gPad.SaveAs(fname)


