#!/usr/bin/env python
#-----------------------------------------------------------------------------
# plotCIspectra.py - example for Roberto
# created 20-Oct-2015 Harrison B. Prosper
#-----------------------------------------------------------------------------
import os, sys
from time import sleep
from histutil import *
from ROOT import *
#-----------------------------------------------------------------------------
LAMBDA  = [10.0, 15.0,      20.0]                # In TeV
COLOR   = [kRed, kOrange+1, kGreen+1]
#-----------------------------------------------------------------------------
def main():
    # load CI library
    gSystem.Load("libCI.so")

    # set up graphics style (from histutil)
    setStyle()

    QCDdir = '../fastNLO/CT10/099'
    CIdir  = '../fastCI/CT10/099'
    
    # -------------------------------------                
    # get list of histograms from qcd.root
    # -------------------------------------            
    rootfile = '%s/qcd.root' % QCDdir
    if not os.path.exists(rootfile):
        hutil.error("plotCIspectra.py",
                    "can't find rootfile %s" % rootfile)
    hfile = TFile(rootfile)
    if not hfile.IsOpen():
        hutil.error("plotCIpectra.py",
                    "can't open rootfile %s" % rootfile)
    histnames = hutil.histogramNames(hfile)
    hfile.Close()

    # create CI spectrum calculator (for first histogram)
    print histnames[0]
    CI  = CISpectrum("CI", "CI", CIdir, histnames[0])

    # --------------------------------------------------------
    filename = 'fig_CI_spectra'
    xmin = 500           # min(pT) GeV
    xmax = 2500          # max(pT) GeV
    ymin =-5.e-4
    ymax = 5.e-5
    
    xx = array('d'); xx.append(xmin); xx.append(xmax)
    yy = array('d'); yy.append(0); yy.append(0)
    horiz = TGraph(2, xx, yy)
    horiz.SetLineWidth(2)
    horiz.SetLineColor(kBlue)
    horiz.GetHistogram().SetAxisRange(xmin, xmax, "X")
    horiz.GetHistogram().SetAxisRange(ymin, ymax, "Y")
    horiz.GetXaxis().SetTitle('Jet p_{T} (GeV)')
    horiz.GetYaxis().SetTitle('CI')

    
    cci = TCanvas(filename,
                  filename,
                  10, 10, 500, 500)
    cci.cd()
    #gPad.SetLogy(kFALSE)
    horiz.Draw('ac')
   
    # Set up a specific CI model
    # A CI model is specified by setting the 6 values of kappa
    # and specifying a value of the mass scale Lambda
    kappa = vector('double')(6, 0)
    kappa[0] = 1
    
    hist = []
    for ii, L in enumerate(LAMBDA):
        l = 1.0/L**2
        color = COLOR[ii]

        # create a histogram containing CI spectrum for
        # specified CI model
        hist.append( CI(l, kappa).Clone() )
        hist[-1].SetAxisRange(xmin, xmax, "X")
        hist[-1].SetLineColor(color)
        hist[-1].SetLineWidth(2)         

        cci.cd()
        hist[-1].Draw('c same')
        cci.Update()
    cci.Update()
    sleep(10)
#-----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print 'ciao!'
