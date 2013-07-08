#include "../interface/FunctionHelpers.h"

// ------------------------------------------------------------------------------------------------
TH1 * integrate1D(TH1 * h, bool normalize) {
	TH1 * ret= (TH1*)h->Clone( Form("%s_cdf", h->GetName() ) );
	ret->SetDirectory(0);
	for(int xx=ret->GetNbinsX()-1; xx>=0; --xx) {
		ret->SetBinContent( xx, ret->GetBinContent(xx) + ret->GetBinContent(xx+1) );
	}
	if( normalize ) { ret->Scale( 1./h->Integral() ); }
	return ret;
}

// ------------------------------------------------------------------------------------------------
TH2 * integrate2D(TH2 * h, bool normalize) {
	TH2 * ret= (TH2*)h->Clone( Form("%s_cdf", h->GetName() ) );
	ret->SetDirectory(0);
	for(int xx=ret->GetNbinsX(); xx>=0; --xx) {
		for(int yy=ret->GetNbinsY()-1; yy>=0; --yy) {
			ret->SetBinContent( xx, yy, ret->GetBinContent(xx,yy) + ret->GetBinContent(xx,yy+1) );
		}
	}
	for(int yy=ret->GetNbinsY(); yy>=0; --yy) {
		for(int xx=ret->GetNbinsX()-1; xx>=0; --xx) {
			ret->SetBinContent( xx, yy, ret->GetBinContent(xx,yy) + ret->GetBinContent(xx+1,yy) );
		}
	}
	if( normalize ) { ret->Scale( 1./h->Integral() ); }
	return ret;
}

// ------------------------------------------------------------------------------------------------
HistoConverter * cdfInv(TH1 * h,double min, double max)
{
	TH1 * hi = integrate1D(h);
	TGraph g(hi);
	TGraph ginv(g.GetN());
	double last = 1.;
	int ilast = 0;
	for(int ip=0; ip<g.GetN(); ++ip) {
		double x = g.GetX()[ip];
		double y = g.GetY()[ip];
		/// std::cout << x << " " << y << std::endl;
		ginv.SetPoint(ip,1.-y,x);
		if( y < last ) { 
			last = y;
			ilast = ip;
			if( ilast == 0 && ip > 1 ) {
				min = x;
			}
		}
	}
	max = g.GetX()[ilast];
	
	HistoConverter * invg = new LinGraphToTF1(Form("%s_inv",hi->GetName()), &ginv, 0., min, 1., max );
	std::cout << min << " " << max << " " << invg->eval(max) << " " << invg->eval(g.GetX()[g.GetN()-1]) << " " << invg->eval(1.) << std::endl;
	delete hi;
	return invg;
}

// ------------------------------------------------------------------------------------------------
HistoConverter * cdf(TH1 * h,double min, double max)
{
	TH1 * hi = integrate1D(h);
	TGraph g(hi);
	TGraph ginv(g.GetN());
	double last = 1.;
	int ilast = 0;
	for(int ip=0; ip<g.GetN(); ++ip) {
		double x = g.GetX()[ip];
		double y = g.GetY()[ip];
		if( y < last ) { 
			last = y;
			ilast = ip;
			if( ilast == 0 && ip > 1 ) {
				min = x;
			}
		}
		ginv.SetPoint(ip,x,1.-y);
	}
	max = g.GetX()[ilast];
	

	HistoConverter * invg = new LinGraphToTF1(Form("%s_dir",hi->GetName()), &ginv, min, 0., max, 1. );
	std::cout << min << " " << max << " " << invg->eval(1.-last) << " " << g.GetY()[g.GetN()-1] << " " << invg->eval(g.GetY()[g.GetN()-1]) << " " << invg->eval(1.) << std::endl;
	delete hi;
	return invg;
}
