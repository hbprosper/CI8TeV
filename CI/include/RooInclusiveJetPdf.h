#ifndef ROOINCLUSIVEJETPDF_H
#define ROOINCLUSIVEJETPDF_H
//---------------------------------------------------------------------------
// Description: model pdf of inclusive jet spectrum averaged
//              over nuisance parameters
// Created: 11-Nov-2014 Harrison B. Prosper
//          09-Feb-2015 HBP add setNumber to allowing picking individual
//                      likelihoods from the ensemble of likelihoods.
//---------------------------------------------------------------------------
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "QCDSpectrum.h"
#include "CISpectrum.h"
#include "Math/Interpolator.h"
//---------------------------------------------------------------------------
class RooInclusiveJetPdf : public RooAbsPdf
{
public:
  RooInclusiveJetPdf() {}
  RooInclusiveJetPdf(const char *name, const char *title,
		     RooArgSet&   _count,
		     RooAbsReal&  _lambda,
		     RooArgSet&   _kappa);

  RooInclusiveJetPdf(const RooInclusiveJetPdf& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const 
  { return new RooInclusiveJetPdf(*this, newname); }
  
  void add(QCDSpectrum& _qcd, CISpectrum& _ci);

  virtual ~RooInclusiveJetPdf() {}

  int getAnalyticalIntegral(RooArgSet& allVars,
  			    RooArgSet& analVars,
  			    const char* /*rangeName*/) const;

  double analyticalIntegral(int code,
  			    const char* /*rangeName*/) const;

  QCDSpectrum* QCD(int c);

  CISpectrum*  CI(int c);
  
  size_t size() { return qcd.size(); }

  void useBinRange(int first=-1, int last=-1);
  void useNumber(int which=0);
  void useAsimov(double lumi, bool yes=true, double l=0);
  
  void initialize();
  
  std::vector<double>& Asimov() { return asimov; }
  std::vector<double>& crossSection(int n);
  std::vector<double>& qcdXsection() { return qcdxsect; }
  		    
 protected:
  double evaluate() const;

  RooListProxy count;
  RooRealProxy lambda;
  RooListProxy kappa;

 private:

  std::vector<QCDSpectrum> qcd;
  std::vector<CISpectrum>  ci;
 
  double smallest;
  double offset;
  int number;
  bool useasimov;
  std::vector<double> asimov;
  std::vector<double> qcdxsect;
  std::vector<double> xsection;
  std::vector<double> data;
  std::vector<double> lngammadata;
  double lngammadatatotal;
  int firstbin;
  int lastbin;
  bool useinterpolation;
  mutable ROOT::Math::Interpolator* interp;

  ClassDef(RooInclusiveJetPdf,1)  
};
//---------------------------------------------------------------------------
#endif
