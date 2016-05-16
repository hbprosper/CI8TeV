#!/usr/bin/env python
#-------------------------------------------------------------------------------
import os, sys, re
from string import *
from glob import glob
from ROOT import *
#-------------------------------------------------------------------------------
HISTNAMES = '''
_0.500_0.500
_0.500_1.000
_1.000_0.500
_1.000_1.000
_1.000_2.000
_2.000_1.000
_2.000_2.000
'''
HISTNAMES = split(strip(HISTNAMES))
MURMUF = [0, 3, 1, 4, 7, 5, 8]
#-------------------------------------------------------------------------------
def xsection(cixsec, l, kappa, which=1):
    bi  = vector('double')(6)
    ai  = vector('double')(6)
    aig = vector('double')(6)

    bij = vector('double')(9)
    aij = vector('double')(9)
    aijg= vector('double')(9)    

    bi4 = vector('double')(4)
    ai4 = vector('double')(4)
    ai4g= vector('double')(4)

    print "aij[3] = ", aij[3]
    
    xsec = [0]*9
    for ii in xrange(9):
        for jj in xrange(6):
            bi[jj] = cixsec.bi(ii, jj, which)
            aig[jj]= cixsec.aig(ii,jj, which)            
            ai[jj] = cixsec.ai(ii, jj, which)

        for jj in xrange(9):
            bij[jj] = cixsec.bij(ii, jj, which)
            aijg[jj]= cixsec.aijg(ii,jj, which)
            aij[jj] = cixsec.aij(ii, jj, which)


        for jj in xrange(4):
            bi4[jj] = cixsec.bi4(ii, jj, which)
            ai4g[jj]= cixsec.ai4g(ii,jj, which)            
            ai4[jj] = cixsec.ai4(ii, jj, which)

        xsec[ii] = CIXsection.evaluate(l, kappa,
                                       bi, aig, ai,
                                       bij,aijg,aij,
                                       bi4,ai4g,ai4)
    return xsec
#-------------------------------------------------------------------------------
def main():
    gSystem.Load('libCI.so')

    Lambda = 10.0
    l = 1.0/Lambda**2
    kappa = vector('double')(6, 0)
    kappa[0] =-1

    histdir = '../CT10/000'
    
    filenames = glob('%s/ci*.txt' % histdir)
    filenames.sort()

    out = open('cixsec.txt', 'w')
    
    for filename in filenames:
        cixsec = CIXsection(filename)

        lo  = xsection(cixsec, l, kappa, 0)
        nlo = xsection(cixsec, l, kappa, 1)

        loxsec = cixsec(l, kappa, 0)
        nloxsec= cixsec(l, kappa, 1)

        record = "==> %s" % filename
        
        out.write('%s\n' % record)     
        print record
        
        record = " %6s %6s\t%12s %12s %10s %10s" % ('mur', 'muf',
                                                 'LO xsec(pb)',
                                                 'NLO xsec(pb)',
                                                 'LO', 'NLO')
        out.write('%s\n' % record)
        print record
        
        muf = 0.50
        mur = 0.50
        for ii in xrange(loxsec.size()):
            if mur > 2:
                muf *= 2
                mur = 0.50
            record = " %6.3f %6.3f\t%12.3e %12.3e %10.3e %10.3e" % (mur, muf,
                                                                 loxsec[ii],
                                                                 nloxsec[ii],
                                                                 lo[ii],
                                                                 nlo[ii])
            out.write('%s\n' % record)
            print record
            
            mur *= 2
    out.close()
#------------------------------------------------------------------------------ 
main()

