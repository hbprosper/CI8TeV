#!/usr/bin/env python
import os, sys, re
from random import normalvariate
from histutil import *
from ROOT import *
#----------------------------------------------------------------------------

xbins = 40
xmin  =-20.0
xmax  = 20.0
 
h1 = mkhist1('h1', 'x', '', xbins, xmin, xmax)
for ii in xrange(1000):
    x = normalvariate(4, 4)
    h1.Fill(x)
    
h2 = mkhist1('h2', 'x', '', xbins, xmin, xmax)
for ii in xrange(1000):
    x = normalvariate(1, 6)
    h2.Fill(x)

h = h1.Clone('h')
h.Reset()
h += 2 * h1
h += 3 * h2
c = TCanvas('c', '', 10, 10, 500, 500)
c.cd()
h.Draw()
c.Update()
gApplication.Run()

    



