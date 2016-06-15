//---------------------------------------------------------------------------
// File: RooInclusiveJetPdf.cc
// Description: compute inclusive jet pdf
// Created: 11-Nov-2014 Harrison B. Prosper
//---------------------------------------------------------------------------
#include <vector>
#include <iostream>
#include <limits>
#include "TMath.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooInclusiveJetPdf.h"
#include "TMath.h"
#include "TRandom3.h"
#include "Math/Random.h"
#include "Math/GSLRndmEngines.h"
//---------------------------------------------------------------------------
ClassImp(RooInclusiveJetPdf)
//---------------------------------------------------------------------------
using namespace std;
namespace {
  ROOT::Math::Random<ROOT::Math::GSLRngMT>* gslrandom = 
    new ROOT::Math::Random<ROOT::Math::GSLRngMT>();
  TRandom3 rand3;
};

RooInclusiveJetPdf::RooInclusiveJetPdf(const char *name, const char *title,
				       RooArgSet&   _count,
				       RooAbsReal&  _lambda,
				       RooArgSet&   _kappa)
  : RooAbsPdf(name, title),
    count("count",   "N",       this, kTRUE, kFALSE), 
    lambda("lambda", "#lambda", this, _lambda),
    kappa("kappa",   "#kappa",  this, kTRUE, kFALSE),
    qcd(vector<QCDSpectrum>()),
    ci(vector<CISpectrum>()),
    smallest(numeric_limits<double>::denorm_min()),
    offset(0),
    number(-1),
    useasimov(false),
    asimov(vector<double>(_count.getSize(), 0)),
    qcdxsect(vector<double>(_count.getSize(), 0)),
    xsection(vector<double>(_count.getSize(), 0)),
    firstbin(0),
    lastbin(_count.getSize()-1),
    useinterpolation(false), // Set false initially to compute likelihood
    interp(0)
{  
  count.add(_count);
  kappa.add(_kappa);
}

RooInclusiveJetPdf::RooInclusiveJetPdf(const RooInclusiveJetPdf& other, 
				       const char* name) 
  : RooAbsPdf(other, name),
    count("count",   this, other.count), 
    lambda("lambda", this, other.lambda),
    kappa("kappa",   this, other.kappa),
    qcd(other.qcd),
    ci(other.ci),
    smallest(other.smallest),
    offset(other.offset),
    number(other.number),
    useasimov(other.useasimov),
    asimov(other.asimov),
    qcdxsect(other.qcdxsect),
    xsection(other.xsection),
    firstbin(other.firstbin),
    lastbin(other.lastbin),
    useinterpolation(other.useinterpolation),
    interp(other.interp)
{}

void RooInclusiveJetPdf::add(QCDSpectrum& _qcd, CISpectrum&  _ci)
{
  qcd.push_back(_qcd);
  ci.push_back(_ci);
}

QCDSpectrum* RooInclusiveJetPdf::QCD(int c)
{
  if ( c < 0 ) return 0;
  if ( c >= (int)qcd.size() ) return 0;
  return &(qcd[c]);
}

CISpectrum* RooInclusiveJetPdf::CI(int c)
{
  if ( c < 0 ) return 0;
  if ( c >= (int)ci.size() ) return 0;
  return &(ci[c]);
}

vector<double>& RooInclusiveJetPdf::crossSection(int c)
{
  for(size_t ii=0; ii < xsection.size(); ii++)  xsection[ii] = 0;
  if ( c < 0 ) return xsection;
  if ( c >= (int)ci.size() ) return xsection;

  // CI parameters
  double l = (double)lambda;
  std::vector<double> k(6);
  for(int ii=0; ii < 6; ii++)
    k[ii] = dynamic_cast<RooRealVar*>(&kappa[ii])->getVal();

  for(size_t ii=0; ii < xsection.size(); ii++)
    {
      qcdxsect[ii] = qcd[c](ii);
      xsection[ii] = qcdxsect[ii] + ci[c](l, k, ii);
    }
  return xsection;
}

void RooInclusiveJetPdf::useNumber(int which)
{
  number = which;
  if ( number > (int)qcd.size()-1 ) number =-1;
}

void RooInclusiveJetPdf::useBinRange(int first, int last)
{
  if ( first < 0 ) first = 0;
  if ( first > count.getSize()-1 ) first = count.getSize()-1;
  if ( last  < 0 ) last  = count.getSize()-1;
  if ( last  > count.getSize()-1 ) last  = count.getSize()-1;
  if ( first > last ) first = last;
  firstbin = first;
  lastbin  = last;
  cout << "RooInclusiveJetPdf: bin range["
       << firstbin << "..." << lastbin << "]"
       << endl;
}

void RooInclusiveJetPdf::useAsimov(double lumi, bool yes, double l)
{
  useasimov = yes;
  if ( ! useasimov ) return;

  std::vector<double> k(6);
  if ( l > 0 )
    for(int c=0; c < 6; c++)
      k[c] = dynamic_cast<RooRealVar*>(&kappa[c])->getVal();
  
  // Note: the first element is presumed to be the
  // nominal cross section.
  double scale = 1000*lumi;
  for(size_t ii=0; ii < asimov.size(); ii++)
    {
      double xsec = qcd[0](ii);
      if ( l > 0 ) xsec += ci[0](l, k, ii);
      asimov[ii] = scale * xsec;
    }
}

void RooInclusiveJetPdf::initialize()
{
  // ----------------------------------------------------
  // compute pdf at different values of Lambda
  // ----------------------------------------------------
  int nbins = lastbin-firstbin+1;
  data.clear();
  lngammadata.clear();
  lngammadatatotal = 0; 
  int jj=0;
  for(int ii=firstbin; ii <= lastbin; ii++)
    {
      if ( useasimov )
	data.push_back(asimov[ii]);
      else
	data.push_back(dynamic_cast<RooRealVar*>(&count[ii])->getVal());
      lngammadatatotal += data[jj];
      lngammadata.push_back(TMath::LnGamma(data[jj]+1));
      jj++;
    }
  lngammadatatotal = TMath::LnGamma(lngammadatatotal+1);
  
  int npts=200;
  vector<double> x;
  vector<double> y;
  RooRealVar* lambdaval = dynamic_cast<RooRealVar*>(lambda.absArg());
  assert(lambdaval);
  double xmin = lambdaval->getMin();
  double xmax = lambdaval->getMax();
  double xstep= (xmax-xmin)/npts;
  
  useinterpolation = false;
  double total = 0;
  for(int c=0; c <= npts; c++)
    {
      double xp = xmin + c * xstep; lambdaval->setVal(xp);
      double yp = evaluate();
      x.push_back(xp);
      y.push_back(yp);
      total += yp;
    }
  useinterpolation = true;
  
  total *= xstep;
  for(size_t c=0; c < y.size(); c++) y[c] /= total;
 
  if ( interp ) delete interp;
  interp = new ROOT::Math::Interpolator(x, y,
					ROOT::Math::Interpolation::kLINEAR);

}

double RooInclusiveJetPdf::evaluate() const
{
  double l = (double)lambda;
  if ( useinterpolation ) return interp->Eval(l);
  
  // CI parameters
  std::vector<double> k(6);
  for(int c=0; c < 6; c++)
    k[c] = dynamic_cast<RooRealVar*>(&kappa[c])->getVal();
  
  // integrated likelihood
  int nbins = lastbin-firstbin+1;
  vector<double> p(nbins,0);
  double y  = 0;
  if ( number > -1 )
    {
      int jj = 0;
      double sum=0;
      for(int ii=firstbin; ii <= lastbin; ii++)
	{
	  p[jj] = qcd[number](ii) + ci[number](l, k, ii);
	  sum += p[jj];
	  jj++;
	}
      double nlnp = lngammadatatotal;
      for(size_t j=0; j < p.size(); j++)
	nlnp += data[j] * log(p[j]/sum) - lngammadata[j];
      y = exp(nlnp);      
    }
  else
    {
      // Note: the first element is presumed to be the
      // nominal cross section so we skip it.
      for(size_t c=1; c < qcd.size(); c++)
	{
	  int jj=0;
	  double sum=0;
	  for(int ii=firstbin; ii <= lastbin; ii++)
	    {
	      p[jj] = qcd[c](ii) + ci[c](l, k, ii);
	      sum += p[jj];
	      jj++;
	    }
	  double nlnp = lngammadatatotal;
	  for(size_t j=0; j < p.size(); j++)
	    nlnp += data[j] * log(p[j]/sum) - lngammadata[j];
	  y += exp(nlnp);
	}
      y /= qcd.size()-1;
    }
  if ( y != y || y <= 0 ) 
    return 0;
  else
    return y;
}

int RooInclusiveJetPdf::getAnalyticalIntegral(RooArgSet& allVars,
					      RooArgSet& analVars, 
					      const char* /*rangeName*/)
  const 
{
  if (matchArgs(allVars, analVars, count)) 
    return 1;
  else
    return 0;
}

double RooInclusiveJetPdf::analyticalIntegral(int code, 
					      const char* /*rangeName*/)
  const 
{
  if(code==1)
    {
      // already normalized!
      return 1;
    } 
  else
    return 0; 
}
