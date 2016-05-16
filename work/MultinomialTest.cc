//--------------------------------------------------------------------
// File: MultinomialTest.cc
// Description: Test multinomial computation in
//				'RooMultinomial.cc' using CI histograms		
// Created: 10-June-2014 Greg Myers
//--------------------------------------------------------------------
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <limits>
#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TMath.h"
#include "RooFit.h"
#include "RooFormulaVar.h"
#include "Math/Random.h"
#include "Math/GSLRndmEngines.h"

using namespace std;

TH1D* MultinomialTest(	TH1F hdata,		TH1F hqcd, 
						TH1F hA,		TH1F hB,
						TH1F hAprime,	TH1F hBprime,
						TH1F hAprimeln,	TH1F hBprimeln )
{
	cout<<"---------------------------------"<<endl;
	cout<<"Testing multinomial calculation: " <<endl;
	cout<<"---------------------------------"<<endl;
	
	double Ntotal		= 0.;//hdata.Integral();
	double sigma_tot	= 0.;
	double l			= 0.01; //TeV^-2
	
	vector<double> N, QCD, A, B, Aprime, Bprime, Aprimeln, Bprimeln;
	vector<double> pT;
	vector<double> sigma;
	double result = 0.0;
		
	//read in histogram values and calculate sigmas
	
	const int nbins = hdata.GetNbinsX();
	double ptbins[nbins];
	for(int i=0;i<nbins;i++)
	{
		pT.push_back( hdata.GetBinLowEdge(i+1) );
		ptbins[i] = hdata.GetBinLowEdge(i+1);
		N.push_back( hdata.GetBinContent(i+1) );
		QCD.push_back( hqcd.GetBinContent(i+1) );
		A.push_back( hA.GetBinContent(i+1) );
		B.push_back( hB.GetBinContent(i+1) );
		Aprime.push_back( hAprime.GetBinContent(i+1) );
		Bprime.push_back( hBprime.GetBinContent(i+1) );
		Aprimeln.push_back( hAprimeln.GetBinContent(i+1) );
		Bprimeln.push_back( hBprimeln.GetBinContent(i+1) );
		
		double c =	QCD[i] 
					+ (B[i] + Bprime[i]*log(1000.0/sqrt(l)) - Bprimeln[i])*l
					+ (A[i] + Aprime[i]*log(1000.0/sqrt(l)) - Aprimeln[i])*l*l;
		
		//cout<<"c = "<<c<<endl;
		Ntotal += N[i];
		sigma_tot += c;
		sigma.push_back( c );
			
	}
	ptbins[20]=2500.0;
	//cout<<"##### Ntotal = " << Ntotal <<" ####"<<endl;
	//cout<<"##### sigma_tot = " << sigma_tot <<" ####"<<endl;
	
//--------------------------------------------------------------------
	TH1D* hist_nlnp = new TH1D( "hist_nlnp", "hist_nlnp",nbins,ptbins);	
	hist_nlnp->GetXaxis()->SetTitle("p_{T} (GeV)");
	hist_nlnp->GetXaxis()->SetRangeUser(507.,2500.);
	hist_nlnp->GetYaxis()->SetTitle("Nlog#left(#frac{#sigma}{#sigma_{tot}}#right)");
	hist_nlnp->GetYaxis()->SetRangeUser(-7.0e5,100.0);
	
//Taken from RooMultinomial.cc:	
	
	// normalize sigmas
	for(int i=0; i < nbins; i++){sigma[i] /= sigma_tot;}
	
	// compute multinomial
	double u = TMath::LnGamma(Ntotal+1);
	cout<< "u = " << u <<endl;
	double y = 0.0;
	double z = 0.0;
	for(int i=0; i < nbins; i++)
    {
		double lng = TMath::LnGamma(N[i]+1);
		y += lng;
		double nlnp = N[i] * log(sigma[i]);
		//double lnp = log(sigma[i]);
		hist_nlnp->SetBinContent(i+1, nlnp);
		z += nlnp;
		cout<< "y = " << y <<setw(15)<< "z = " << z <<endl;
    }
	result = z + u - y;
	if ( result != result){exit(0);}
	
	cout << "log(result) = " << result <<endl;
	result = exp(result);
	
	cout << "result = "<< result <<endl;
	
	cout<<"---------------------------------"<<endl;
	cout<<"     end of Multinomial Test     "<<endl;
	cout<<"---------------------------------"<<endl;
	
	
	return hist_nlnp;
}

