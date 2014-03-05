#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-d","--dir")
parser.add_option("-s","--spinCats",type="int",default=5)
parser.add_option("-y","--yaxis",default="-1.5,3.5")
parser.add_option("-u","--unblind",default=False,action="store_true")
parser.add_option("-b","--batch",default=False,action="store_true")
parser.add_option("-o","--outfile",default="chcomp")
parser.add_option("-S","--sqrtS",default="comb")
parser.add_option("-m","--mass",type="float",default=120)
(options,args) = parser.parse_args()

import ROOT as r
import sys,os

if options.batch:
	r.gROOT.SetBatch()

def findQuantile(pts,cl):

	#gr is a list of r,nll
	# start by walking along the variable and check if crosses a CL point
	crossbound = [ pt[1]<=cl for pt in pts ]
	rcrossbound = crossbound[:]
	rcrossbound.reverse()

	minci = 0
	maxci = len(crossbound)-1
	min = pts[0][0]
	max = pts[maxci][0]

	for c_i,c in enumerate(crossbound): 
		if c : 
			minci=c_i
			break
	
	for c_i,c in enumerate(rcrossbound): 
		if c : 
			maxci=len(rcrossbound)-c_i-1
			break

	if minci>0: 
		y0,x0 = pts[minci-1][0],pts[minci-1][1]
		y1,x1 = pts[minci][0],pts[minci][1]
		min = y0+((cl-x0)*y1 - (cl-x0)*y0)/(x1-x0)
		
	if maxci<len(crossbound)-1: 
		y0,x0 = pts[maxci][0],pts[maxci][1]
		y1,x1 = pts[maxci+1][0],pts[maxci+1][1]
		max = y0+((cl-x0)*y1 - (cl-x0)*y0)/(x1-x0)

	return min,max

def plot1DNLL(file,name=""):

	x = 'r'

	clean_graphs=[]

	tf = r.TFile(file)
	tree = tf.Get('limit')
	gr = r.TGraph()
	gr.SetName('gr_%s_%s'%(x,name))

	res=[]
	for i in range(tree.GetEntries()):
		tree.GetEntry(i)
		xv = getattr(tree,x)
		#if tree.quantileExpected==1: continue
		if tree.deltaNLL<0: continue
		if xv in [re[0] for re in res]: continue
		if tree.deltaNLL>100 or tree.deltaNLL<-100: continue
		if tree.deltaNLL != tree.deltaNLL: continue
		if xv<-1 and tree.deltaNLL==0: continue
		res.append([xv,2.*tree.deltaNLL])
	res.sort()
	minNLL = min([re[1] for re in res])
	for re in res: 
		#if options.correctNLL and re[1]==0.: re[1]=-1
		re[1]-=minNLL
	
	p=0
	for re, nll in res: 
		if nll>=0.:
			gr.SetPoint(p,re,nll)
			p+=1

	m,m1 = findQuantile(res,0);
	#m = 1.
	#m1 = 1.
	l,h  = findQuantile(res,1);
	l2,h2  = findQuantile(res,4);

	"""
	if options.shiftNLL:
		shiftNLL(gr,m)
		shift = options.shiftNLL - m
		m += shift
		m1 += shift
		l += shift
		l2 += shift
		h += shift
		h2 += shift

	if options.smoothNLL:
		smoothNLL(gr,res)
		clean_graphs.append(gr)
	else:
	"""
	clean_graphs.append(gr)

	xmin = m
	eplus = h-m
	eminus = m-l
	eplus2 = h2-m
	eminus2 = m-l2

	print "%15s - %4.2f +%4.2g -%4.2g" % ( name, xmin, eplus , eminus )
	#print "%s : %1.4f +%1.3g -%1.3g (2sig) +%1.3g -%1.3g" % ( options.names[k], xmin, eplus , eminus, eplus2, eminus2 )

	#canvTemp = r.TCanvas()
	#gr.Draw("ALP")
	#raw_input("?\n")

	return [xmin,eplus,eminus,eplus2,eminus2]

files={}
files['Data'] = []
files['SM'] = []
files['GravGG'] = []
files['Grav50GG50QQ'] = []
files['GravQQ'] = []
scale={}
scale['Data'] = 1.
scale['SM'] = 1.
scale['GravGG'] = 1.4
scale['Grav50GG50QQ'] = 1.265
scale['GravQQ'] = 1.019
graphs={}
opts={}
opts['Data'] = 						[r.kBlack,	r.kFullCircle,			"Observed"]
opts['SM'] = 							[r.kRed,		r.kFullSquare,			"X#rightarrow#gamma#gamma 0^{+}"]
opts['GravGG'] = 					[r.kBlue,		r.kFullTriangleUp,	"X#rightarrow#gamma#gamma 2^{+}_{m}(100%gg)"]
opts['Grav50GG50QQ'] = 		[r.kMagenta,r.kFullDiamond,			"X#rightarrow#gamma#gamma 2^{+}_{m}(50%gg,50%qq)"]
opts['GravQQ'] = 					[r.kGreen,	r.kFullTriangleDown,"X#rightarrow#gamma#gamma 2^{+}_{m}(100%qq)"]
dummyHist = r.TH2F('dummy','',1,0.,1.,1,float(options.yaxis.split(',')[0]),float(options.yaxis.split(',')[1]))
dummyHist.GetXaxis().SetTitle('|cos#theta*|')
dummyHist.GetYaxis().SetTitle('#sigma/#sigma_{SM}')
x=[0.1,0.2875,0.4625,0.65,0.875]

canv = r.TCanvas()
canv.SetGrid()

dummyHist.Draw("AXISG")
dummyHist.SetStats(0)
leg = r.TLegend(0.11,0.65,0.4,0.89)
leg.SetFillColor(0)
leg.SetLineColor(0)

fkeys = files.keys()
fkeys.sort()
print fkeys

for key in fkeys:
	graphs[key] = r.TGraphAsymmErrors()
	graphs[key].SetName(key)
	graphs[key].SetLineWidth(2)
	graphs[key].SetLineColor(opts[key][0])
	graphs[key].SetMarkerColor(opts[key][0])
	graphs[key].SetMarkerStyle(opts[key][1])
	if key=='Data': leg.AddEntry(graphs[key],opts[key][2],"EP")
	else: leg.AddEntry(graphs[key],opts[key][2],"LP")
	for p,cat in enumerate(range(options.spinCats)):
		temp_files = []
		temp_files.append(('%s/higgsCombineSMto%sCat%d.MultiDimFit.mH%6.2f.root'%(options.dir,key,cat,options.mass)).replace('.00',''))
		res = plot1DNLL(temp_files[0],'%s_cat%d'%(key,cat))
		#graphs[key].SetPoint(p,x[p],scale[key]*res[0])
		graphs[key].SetPoint(p,x[p],res[0])
		if key=='Data': graphs[key].SetPointError(p,0.,0.,res[2],res[1])
		else: graphs[key].SetPointError(p,0.,0.,0.,0.)
		#array.append(r.TFile('%s/higgsCombineSMto%sCat%d.MultiDimFit.mH120.root'%(options.dir,key,cat)))
	
	if key=='Data': 
		if options.unblind: 
			graphs[key].Draw("PEsame")	
	else: 
		graphs[key].Draw("LPsame")

if options.unblind: graphs['Data'].Draw("PEsame")
leg.Draw("same")
pt = r.TPaveText(0.1,0.91,0.45,0.99,"NDC");
pt.SetTextAlign(12);
pt.SetTextSize(0.04);
pt.SetFillColor(0);
pt.AddText("CMS Preliminary");
pt.SetBorderSize(0);
pt2 = r.TPaveText(0.75,0.90,0.9,0.99,"NDC");
pt2.SetTextAlign(32);
pt2.SetTextSize(0.04);
pt2.SetFillColor(0);
if options.sqrtS=="comb":
	pt2.AddText("#splitline{#sqrt{s} = 7 TeV, L = 5.1 fb^{-1}}{#sqrt{s} = 8 TeV, L = 19.7 fb^{-1}}");
elif options.sqrtS=="7" or options.sqrtS=="7TeV":
	pt2.AddText(" #sqrt{s} = 7 TeV, L = 5.1 fb^{-1}");
elif options.sqrtS=="8" or options.sqrtS=="8TeV":
	pt2.AddText(" #sqrt{s} = 8 TeV, L = 19.7 fb^{-1}");
else:
	sys.exit("Invalid sqrtS"+options.sqrtS)
pt2.SetBorderSize(0);
pt.Draw("same")
pt2.Draw("same")
canv.Update()
canv.Modified()
if not options.batch:
	raw_input('Good?\n')
canv.Print('%s.pdf'%options.outfile)
canv.Print('%s.png'%options.outfile)
