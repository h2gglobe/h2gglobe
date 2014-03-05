import ROOT as r

r.gROOT.ProcessLine(".L $CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so")

f7 = r.TFile("CMS-HGG_spin_7TeV_data_collapsed.root")
f8 = r.TFile("CMS-HGG_spin_8TeV_data_collapsed.root")
fP = r.TFile("unblinding_v1/hgg/hgg.inputsig_8TeV_spin.root")


for f in [f7,f8]:
	dCount=0.
	smCount = 0. 
	ggCount = 0. 
	qqCount = 0.
	print f.GetName()
	ws = f.Get('cms_hgg_workspace')
	for c in range(4):
		dEnts=0.
		dEnts = ws.data('data_mass_cat%d'%c).sumEntries()
		dCount += dEnts 
		smEnts = 0.
		smArr = ['ggh','vbf','wh','zh','tth']
		if f==fP: smArr = ['ggh','vbf','wzh','tth']
		for p in smArr:
			smEnts += ws.data('sig_%s_mass_m125_cat%d'%(p,c)).sumEntries()
		smCount += smEnts
		ggEnts = ws.data('sig_gg_grav_mass_m125_cat%d'%c).sumEntries()
		ggCount += ggEnts
		qqEnts = ws.data('sig_qq_grav_mass_m125_cat%d'%c).sumEntries()
		qqCount += qqEnts

		print 'cat', c, ' -- D:', dEnts, 'SM:', smEnts, 'GG:', ggEnts, 'QQ:', qqEnts
	print 'TOTAL: -- D:', dCount, 'SM:', smCount, 'GG:', ggCount, 'QQ:', qqCount




