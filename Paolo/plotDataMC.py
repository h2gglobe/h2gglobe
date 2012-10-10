#! /usr/bin/env python

from optparse import OptionParser, make_option
import sys
import os
import math
import time
import subprocess as sub
from array import array
from ROOT import *
#print "done."

def main(options,args):
    # if you want to run in batch mode
    ROOT.gROOT.SetBatch()

    plots = [ "pho1_pt", "pho1_eta", "pho2_pt", "pho2_eta", "mass", "all_mass", "nvtx" ]
    datasets = {}
    datasets["0_data"] = [ "Data" ]
    datasets["1_prompt_prompt"] = [ "diphojet_8TeV", "dipho_Box_25_8TeV", "dipho_Box_250_8TeV" ]
    datasets["2_prompt_fake"] = [ "gjet_20_8TeV_pf", "gjet_40_8TeV_pf" ]
    datasets["3_fake_fake"] = [ "qcd_30_8TeV_pf", "qcd_40_8TeV_pf" ]
    datasets["4_DY"] = [ "DYJetsToLL" ]
    datasets["5_signal"] = [ "ggh_m125_8TeV", "vbf_m125_8TeV", "wzh_m125_8TeV", "tth_m125_8TeV" ]

    color = {}
    color["0_data"] = kBlack
    color["1_prompt_prompt"] = kBlue
    color["2_prompt_fake"] = kGreen
    color["3_fake_fake"] = kRed
    color["4_DY"] = kCyan
    color["5_signal"] = kYellow

    input = TFile( options.infile )
    input.Print()
    ROOT.gStyle.SetOptFit(11111)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    
    if os.path.isdir(options.outdir) is False:
        os.mkdir(options.outdir)

    for plot in plots:
        for cat in range(0,10):
            histos = {}
            mc=THStack("hs","")
            a = TLegend(0.7,0.7,0.98,0.98)
            a.SetFillColor(0)
            a.SetBorderSize(0)
            for type in sorted(datasets.keys()):
                isample=0
                h =  input.Get ( "%s_cat%d_%s" % (plot,cat,datasets[type][0]) )
                histos["%s_cat%d_%s" % (plot,cat,type) ] = h.Clone( "%s_cat%d_%s" % (plot,cat,type) )
                histos["%s_cat%d_%s" % (plot,cat,type) ].SetFillColor(color[type])
                histos["%s_cat%d_%s" % (plot,cat,type) ].SetLineColor(color[type])
                histos["%s_cat%d_%s" % (plot,cat,type) ].SetMarkerColor(color[type])
                #                histos["%s_cat%d_%s" % (plot,cat,type) ].Sumw2()
                for dataset in datasets[type][1:]:
                    h1 =  input.Get ( "%s_cat%d_%s" % (plot,cat,dataset) )
                    histos["%s_cat%d_%s" % (plot,cat,type) ].Add(h1)
                if (type != "0_data" and type != "5_signal"):
                    mc.Add(histos["%s_cat%d_%s" % (plot,cat,type) ])
                    a.AddEntry(histos["%s_cat%d_%s" % (plot,cat,type)],"%s=%7.1f" % (type,histos["%s_cat%d_%s" % (plot,cat,type)].Integral()),"F")
                elif (type == "0_data") :
                    a.AddEntry(histos["%s_cat%d_%s" % (plot,cat,type)],"%s=%7.1f" % (type,histos["%s_cat%d_%s" % (plot,cat,type)].Integral()),"PL")
                elif (type == "5_signal") :
                    a.AddEntry(histos["%s_cat%d_%s" % (plot,cat,type)],"%s=%7.1f" % (type,histos["%s_cat%d_%s" % (plot,cat,type)].Integral()),"F")

            canv = TCanvas("%s_cat%d" % (plot,cat), "%s_cat%d" % (plot,cat), 800,600)
            histos["%s_cat%d_%s" % (plot,cat,"0_data")].SetMarkerStyle(20)
            histos["%s_cat%d_%s" % (plot,cat,"0_data")].SetMarkerSize(1.3)
            histos["%s_cat%d_%s" % (plot,cat,"0_data")].GetXaxis().SetTitle(plot)
            histos["%s_cat%d_%s" % (plot,cat,"0_data")].Draw("PE")
            mc.Draw("HSAME")
            histos["%s_cat%d_%s" % (plot,cat,"5_signal")].Draw("HSAME")
            histos["%s_cat%d_%s" % (plot,cat,"0_data")].Draw("PESAME")
            a.Draw()
            canv.SaveAs( os.path.join(options.outdir,"%s_cat%d.%s" % (plot,cat,"png") ) )


if __name__ == "__main__":
    parser = OptionParser(option_list=[
        make_option("-i", "--infile",
                    action="store", type="string", dest="infile",
                    default="histograms_dataMC.root",
                    help="", metavar=""
                    ),
        make_option("-d", "--outdir",
                    action="store", type="string", dest="outdir",
                    default="plots",
                    help="", metavar=""
                    )
        ])
    
    (options, args) = parser.parse_args()
    
    main( options, args) 
