#!/usr/bin/env python

import os
import csv

import PlotUtils

import sys
#----------------------------------------------------------------------
# parameters
#----------------------------------------------------------------------

# scaling factors for plotting the SM cross section
# (just for illustration purposes)
smCrossSectionScalingFactors = [ 1, 3, 5, 10 ]

# phobicCrossSectionScalingFactors = [1,5]
phobicCrossSectionScalingFactors = smCrossSectionScalingFactors[:]

classicMode = False

#----------------------------------------------------------------------

def makeGraphLabelOnRight(graph, minMass, maxMass, text):

    lastY = PlotUtils.getGraphYvalues(graph)[-1]
    
    label = ROOT.TLatex(maxMass + 0.01 * (maxMass - minMass),
                        lastY,
                        text)
    
    label.SetTextAlign(12)

    return label

#----------------------------------------------------------------------

def makePlot(csvFile, relative, includeExpected = True, fermiophobic = None, ymax = None,
             ignoreMassesExpected = None,
             exclude_null_expected = False,
             ):
    """ @param relative if true, the exclusion of the ratios
        (relative to the inputs given) are plotted. If False,
        these ratios are multiplied by the SM cross sections
        and branching ratios into gamma gamma 
    """
            
    masses, observedValues, expectedValues, \
           expected_minus_2_sigma_values, \
           expected_minus_1_sigma_values, \
           expected_plus_1_sigma_values, \
           expected_plus_2_sigma_values = PlotUtils.readCSV(csvFile, includeExpected)

    massesExpected = masses[:]

    # compile an overall list of mass points to be excluded for the expected limits
    expectedMassExclusionList = []

    # mass points to be excluded explicitly specified
    if ignoreMassesExpected != None:
        expectedMassExclusionList.extend(ignoreMassesExpected)

    # auto-exclusion of those which are zero
    if exclude_null_expected:
        # find those mass points for which the expected limit is zero
        for mass, expectedLimit in zip(massesExpected, expectedValues):
            # TODO: maybe add some tolerance
            if expectedLimit == 0:
                expectedMassExclusionList.append(mass)

    if expectedMassExclusionList:    
        # at least one mass point to be excluded (for the expected
        # limits)

        tolerance = 1e-4

        for mass in expectedMassExclusionList:
            for i, value in enumerate(massesExpected):
                # print "TESTING MASS",value
                if abs(mass - value) < tolerance:
                    # print "EXCLUDING",value
                    del massesExpected[i]
                    del expectedValues[i]
                    del expected_minus_2_sigma_values[i]
                    del expected_minus_1_sigma_values[i]
                    del expected_plus_1_sigma_values[i]
                    del expected_plus_2_sigma_values[i]

                    break

    if not relative:

        # if we're plotting the absolute cross sections, we
        # need to know whether this is Standard Model or Fermiophobic
        assert(fermiophobic != None)

        if fermiophobic:
            typeName = "FP"
        else:
            typeName = "SM"

        observedValues = PlotUtils.multiplyArrayByXsectAndBR(masses, observedValues, fermiophobic)
        expectedValues = PlotUtils.multiplyArrayByXsectAndBR(massesExpected, expectedValues, fermiophobic)

        for values in (expected_minus_2_sigma_values,
                    expected_minus_1_sigma_values,
                    expected_plus_1_sigma_values,
                    expected_plus_2_sigma_values):
            if values == None:
                continue

            # print locals()[name][-1]
            newValues = PlotUtils.multiplyArrayByXsectAndBR(massesExpected, values, fermiophobic)
            # print locals()[name][-1]

            for i,j in enumerate(newValues):
                values[i] = j

        # print "expected_plus_2_sigma_values",expected_plus_2_sigma_values

        if options.printAbsExclusion:
            print "limits on cross sections times branching ratio [pb]:"
            print "  mass        exp     obs"
            for mass,exp,obs in zip(masses, expectedValues, observedValues):
                print "  %.1f GeV: %7.4f %7.4f" % (mass, exp, obs)

    #----------------------------------------
    # produce the plot
    #----------------------------------------

    grObserved = PlotUtils.makeGraphFromArrays(masses, observedValues)
    gcSaver.append(grObserved)

    # make the cross section graphs for the fixed scaling values
    # (only if we're plotting the absolute cross section
    # exclusions)

    xsectGraphs = []

    if fermiophobic:
        crossSectionScalingFactors = phobicCrossSectionScalingFactors
    else:
        crossSectionScalingFactors = smCrossSectionScalingFactors


    
    if not relative:
        # get interpolated cross sections at masses
        # where we calculated the exclusion
        smXsect = [ PlotUtils.getXsectTimesBR(mass, fermiophobic) for mass in masses ]

        if options.printXsect:
            print "used cross sections times branching ratio:"
            for mass,xsectTimesBR in zip(masses,smXsect):
                print "  %.1f GeV: " % mass,xsectTimesBR

        # color = ROOT.TColor(38, )
        color = ROOT.gROOT.GetColor(38)
        color.SetRGB(0, 0.6, 1)

        # create but do not draw yet the cross section graphs 

        for scalingFactor in crossSectionScalingFactors:
            gr = PlotUtils.makeGraphFromArrays(masses, [ xsect * scalingFactor for xsect in smXsect])
            gcSaver.append(gr)
            
            xsectGraphs.append(gr)

            gr.SetLineWidth(4)
            gr.SetLineStyle(ROOT.kDashed)
            # gr.SetLineColor(ROOT.kBlue)
            
            
            gr.SetLineColor(38)
    #----------------------------------------
    # determine the y scale
    #----------------------------------------

    if ymax == None:

        ymax = max(observedValues) 

        if includeExpected:
            ymax = max(ymax, max(expectedValues))

        if xsectGraphs:
            ymax = max(ymax, max([ PlotUtils.getGraphMaximumY(gr) for gr in xsectGraphs]))

        ymax *= 1.1

        if not relative:
            if fermiophobic:
                # fix the y scale by hand in order not to
                # stretch it too much because of large
                # scaling factors for the theoretical expectation
                ymax = 0.1
            

    #----------------------------------------

    # create a dummy histogram to set the x range
    hdummy = ROOT.TH1F("hdummy","",1,min(masses),max(masses))
    gcSaver.append(hdummy)

    hdummy.SetMaximum(ymax)
    hdummy.Draw("AXIS")

    ROOT.gStyle.SetOptTitle(0)

    #----------------------------------------
    # sigma bands
    #----------------------------------------
    if not classicMode and includeExpected and \
       expected_minus_2_sigma_values != None and \
       expected_minus_1_sigma_values != None and \
       expected_plus_1_sigma_values != None and \
       expected_plus_2_sigma_values != None:

        for sigmas, color in (
            (2,ROOT.kYellow),
            (1,ROOT.kGreen)):

            plusValues = locals()['expected_plus_%d_sigma_values' % sigmas]
            minusValues = locals()['expected_minus_%d_sigma_values' % sigmas]

            # TODO: this (making the filled graph) should be made a function in PlotUtils
            # (is also used in ~/2011-09-vbf-cl/mystuff/random-experiments/plot.py)
            xvaluesForFilledGraph = massesExpected[:]
            yvaluesForFilledGraph = plusValues[:]

            xvaluesForFilledGraph.extend(massesExpected[::-1])
            yvaluesForFilledGraph.extend(minusValues[::-1])

            # just to be sure: duplicate first point
            xvaluesForFilledGraph.append(xvaluesForFilledGraph[0])
            yvaluesForFilledGraph.append(yvaluesForFilledGraph[0])

            grFilled = PlotUtils.makeGraphFromArrays(xvaluesForFilledGraph, yvaluesForFilledGraph)

            gcSaver.append(grFilled)

            grFilled.SetFillColor(color)
            grFilled.Draw("F")


    if classicMode:

        #----------------------------------------
        # draw a filled area line for the observed limit
        #----------------------------------------
        # add points at the corners of the plot
        # and add the first point at the end
        # to get a filled polygon, see http://root.cern.ch/root/roottalk/roottalk00/0256.html

        xvaluesForFilledGraph = masses[:]
        yvaluesForFilledGraph = observedValues[:]

        xvaluesForFilledGraph.append(masses[-1]) ; yvaluesForFilledGraph.append(ymax)
        xvaluesForFilledGraph.append(masses[0])  ; yvaluesForFilledGraph.append(ymax)

        xvaluesForFilledGraph.append(masses[0])  ; yvaluesForFilledGraph.append(observedValues[0])

        grFilled = PlotUtils.makeGraphFromArrays(xvaluesForFilledGraph, yvaluesForFilledGraph)

        gcSaver.append(grFilled)

        grFilled.SetFillColor(ROOT.kOrange)
        grFilled.Draw("F")

    else:
        # draw a plain line for the observed
        grObserved = PlotUtils.makeGraphFromArrays(masses, observedValues)
        gcSaver.append(grObserved)
        grObserved.SetLineWidth(4)
        grObserved.Draw("C")
        #grObserved.Draw("L")

    #----------------------------------------
    # expected graph
    #----------------------------------------

    if includeExpected:
        grExpected = PlotUtils.makeGraphFromArrays(massesExpected, expectedValues)
        gcSaver.append(grExpected)
        
        grExpected.SetLineStyle(ROOT.kDashed)
        grExpected.SetLineWidth(4)
        grExpected.SetLineColor(ROOT.kRed)
        
        grExpected.Draw("C,same")
        #grExpected.Draw("L,same")

        label = makeGraphLabelOnRight(grExpected, min(masses), max(masses), "BG exp.")
        
        label.SetTextSize(label.GetTextSize() * 0.7)    
        label.Draw()
        gcSaver.append(label)

        
    #----------------------------------------
    # cross section graphs
    #----------------------------------------

    for gr, scalingFactor in zip(xsectGraphs, crossSectionScalingFactors):
        #gr.Draw("L,same")
        gr.Draw("C,same")

        # label the cross section graph
        # get the y position at the end of the graph

        label = makeGraphLabelOnRight(gr, min(masses), max(masses), "%.0f * %s" % (scalingFactor, typeName))
        
        label.SetTextSize(label.GetTextSize() * 0.7)    
        label.Draw()
        gcSaver.append(label)
    


    #----------------------------------------

    ROOT.gStyle.SetOptStat(0)
    hdummy.SetXTitle("m_{H} [GeV/c^{2}]")

    if relative:
        hdummy.SetYTitle("#sigma/#sigma(theo)")
    else:
        hdummy.SetYTitle("#sigma(%s) #times BR(#gamma#gamma) [pb]" % typeName)
        
    ROOT.gPad.SetGrid()

    label = ROOT.TLatex(0.5,0.85,"Excluded at 95% CL.")
    gcSaver.append(label)
    label.SetNDC(1)

    label.SetTextAlign(21)
    label.Draw()

    # make grid lines appear above graph
    hdummy.Draw("AXIGSAME")
    ROOT.gPad.Modified()

                    

#----------------------------------------------------------------------
# main
#----------------------------------------------------------------------
from optparse import OptionParser

parser = OptionParser("""

%prog [options] csv_file

    produces a cross section exclusion plot based on
    the excluded scaling factors of the signal given
    to the confidence level calculations.
    
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
    "--ymax",
    type=float,
    help="manually sepcify the y scale",
    default=None,
    )

parser.add_option(
    "--printXsect",
    action = "store_true",
    help="print the cross section times branching ratio used to convert from ratios to absolute cross sections excluded (for debugging purposes)",
    default = False,
    )


parser.add_option(
    "--printAbsExclusion",
    action = "store_true",
    help="print the excluded cross sections (observed and expected), mainly useful for debugging/comparison purposes",
    default = False,
    )

parser.add_option(
    "--plotPrefix",
    help="prefix for the output eps and png files (default is the empty string) where xsect-{abs,ratio}-exclusion.{png,eps} is appended",
    default="",
    )

parser.add_option(
    "--ignoreMassesExpected",
    help="a comma separated list of masses which should not be included in the plot of the expected limits",
    default=None,
    type = str,
    )

parser.add_option(
    "--exclude-null-expected",
    help="removes those mass points for the expected limit graphs where the expected limit is zero",
    default = False,
    dest="exclude_null_expected",
    action = "store_true", 
    )




(options, ARGV) = parser.parse_args()

#----------------------------------------
if len(ARGV) != 1:
    PlotUtils.optError("expected exactly one non-option argument.")

#----------------------------------------

if options.ignoreMassesExpected != None:
    options.ignoreMassesExpected = options.ignoreMassesExpected.split(',')

    options.ignoreMassesExpected = [ float(x) for x in options.ignoreMassesExpected ]


#----------------------------------------------------------------------

import ROOT
ROOT.gROOT.SetStyle("Plain")
gcSaver = []

makePlot(open(ARGV[0]), relative = True, ymax = options.ymax,
         ignoreMassesExpected = options.ignoreMassesExpected,
         exclude_null_expected = options.exclude_null_expected,

         )
ROOT.gPad.SaveAs(options.plotPrefix + "xsect-ratio-exclusion.eps")
ROOT.gPad.SaveAs(options.plotPrefix + "xsect-ratio-exclusion.png")

makePlot(open(ARGV[0]), relative = False, fermiophobic = options.fermiophobicMode, ymax = options.ymax,
         ignoreMassesExpected = options.ignoreMassesExpected,
         exclude_null_expected = options.exclude_null_expected,
         )

ROOT.gPad.SaveAs(options.plotPrefix + "xsect-abs-exclusion.eps")
ROOT.gPad.SaveAs(options.plotPrefix + "xsect-abs-exclusion.png")


