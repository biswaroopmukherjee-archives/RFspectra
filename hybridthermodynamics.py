# for analysis
import pprint, pickle
import hybrid as hb
import bec1db as db
import pandas as pd
import therpy as tp
import numpy as np
import matplotlib.pyplot as plt
import csv
import scipy
import os
import numpy.ma as ma
cst = tp.cst(sigmaf=0.5)
import urllib.request


class Thermodynamics:
    def __init__(self, z, n, bins=70, tfit_U_range=[500,3000]):
        self.z = z
        self.n = n
        self.bins = bins
        self.tfit_U_range = tfit_U_range
        self.getEFbin()
        self.getPtilde()
        self.getTtilde()

    def getEFbin(self):
        self.EF = cst.n2EF(self.n)
        self.TF = self.EF/cst.kB
        self.U = z2u(self.z)
        binlist = binxy(self.U, self.EF, bins=self.bins)
        self.Ubin, self.EFbin = binlist[0], binlist[1]
        self.nbin = cst.EF2n(self.EFbin,neg=True)
        self.TFbin = self.EFbin/cst.kB

    def getPtilde(self):
        self.P_total = scipy.integrate.trapz(y=self.nbin, x=self.Ubin)
        self.cumulative_integral = scipy.integrate.cumtrapz(y=self.nbin, x=self.Ubin, initial=0)
        self.Pressure = self.cumulative_integral-self.P_total
        self.Ptilde = -self.Pressure/((2/5)*self.nbin*self.EFbin)
        self.Ptildefull = np.interp(x=self.EF, xp=self.EFbin, fp=self.Ptilde)

    def getTtilde(self):
        EOS = MarkEOS()
        maskUbin = (self.Ubin/cst.h>self.tfit_U_range[0]) & (self.Ubin/cst.h<self.tfit_U_range[1])
        self.maskUbin = maskUbin
        maskPtilde = ma.masked_array(self.Ptilde, ~maskUbin)
        maskTF = ma.masked_array(self.TFbin, ~maskUbin)
        self.maskedUbin = ma.masked_array(self.Ubin/cst.h, ~maskUbin)

        maskPtilde_data = maskPtilde.compressed()
        maskTF_data = maskTF.compressed()
        maskTtilde_data = np.interp(x=maskPtilde_data, xp=EOS.Ptilde, fp = EOS.Ttilde)
        self.T_unaveraged = maskTtilde_data* maskTF_data
        self.T = np.nanmean(maskTtilde_data* maskTF_data)
        self.Ttilde = self.T/self.TF
        self.Ttilde[self.Ttilde>5]=np.nan

    def plotTtilde(self, ylims=None):
        plt.plot(1e6 *self.z, self.Ttilde,'.',label='data')
        plt.plot([-300, 300], [0.17,0.17],'r',label=r'T_c')
        plt.legend(loc=9)
        if ylims==None:
            plt.ylim([0,1])
        else:
            plt.ylim(ylims)
        plt.xlim([-150, 150])
        plt.xlabel(r'$z [\mu m]$',fontsize=14)
        plt.ylabel(r'$\widetilde{T}$',fontsize=14)
        plt.title(r'$T_{min} = $'+ '%0.2f' % np.nanmin(self.Ttilde))
        plt.savefig('figures/temp_single', dpi=300)
        plt.show()

    def plotPtilde(self, xlims=None):
        eos = MarkEOS()
        plt.plot(eos.Ttilde, eos.Ptilde,'.',label='MarkEOS')
        plt.plot(self.T/self.TFbin, self.Ptilde,'.',label='single shot data')
        plt.legend(loc=4)
        if xlims==None:
            plt.xlim([0,1])
        else:
            plt.xlim(xlims)
        plt.ylim([0,2.5])
        plt.xlabel(r'$\widetilde{T}=T/T_F$',fontsize=14)
        plt.ylabel(r'$\widetilde{P} = P/(0.4nE_F)$',fontsize=14)
        plt.savefig('figures/markeosfit', dpi=300)
        plt.show()

    def plotUnaveragedT(self, xlims=None):
        plt.plot(self.maskedUbin.compressed(), 1e9*self.T_unaveraged, 'o')
        if xlims==None:
            plt.xlim([0, 6e3])
        else:
            plt.xlim(xlims)
        plt.ylabel('T [nK]')
        plt.xlabel('U [Hz]')
        plt.show()


class MarkEOS():
    def __init__(self):
        if not os.path.exists('PT.csv'):
            urllib.request.urlretrieve('https://www.dropbox.com/s/gjhpemvc3he8nar/PT.csv?dl=1', 'PT.csv')

        self.Ttilde = pd.read_csv('PT.csv', header=None)[0].tolist()
        self.Ptilde = pd.read_csv('PT.csv', header=None)[1].tolist()

    def plotPtilde(self):
        plt.plot(self.Ttilde,self.Ptilde,'.')
        plt.show()



## Helper functions
def z2u(z,f=23.9):
    return (cst.mass*((cst.twopi*f)**2)) *(z**2)/2


def binxy(x,y,bins=200):
    hist, xbin_edges = np.histogram(x,bins=bins)
    positions = np.digitize(x,xbin_edges)
    ybins = np.zeros(len(xbin_edges))
    yerrs = np.zeros(len(xbin_edges))

    for bin_index in np.unique(positions):
        accumulate = np.array([])
        for i in range(len(positions)):
            if positions[i]-1 == bin_index-1:
                accumulate = np.append(accumulate, y[i])
        ybins[bin_index-1] = np.nan_to_num(np.nanmean(accumulate))
        yerrs[bin_index-1] = np.nan_to_num(np.nanstd(accumulate)/np.sqrt(len(accumulate)))

    return xbin_edges, ybins, yerrs
