# for analysis
import pprint, pickle
import pandas as pd
import therpy as tp
import numpy as np
import matplotlib.pyplot as plt
import csv
import scipy
import numpy.ma as ma
cst = tp.cst(sigmaf=0.5)
from tqdm import tqdm_notebook as tqdm
import os
from scipy.optimize import brentq

import hybrid as hb
import bec1db as db
import hybridthermodynamics as hbthermo
import virial


# Helper functions - put in bec1db later on
def save_imagenames(filename,images):
    np.savetxt(filename, images, delimiter=",", fmt='%s')

def read_imagenames(filename):
    return [str(name) for name in np.genfromtxt(filename, dtype='str')]

# For plotting
def custom_pad_khz(array):
    pads = array[:,0:1]-1
    return np.concatenate((pads,array),axis=1)

def custom_pad_ef(array):
    pads = array[:,0:1]-0.1
    return np.concatenate((pads,array),axis=1)

def custom_pad_z(array):
    pads = array[:,0:1]
    return np.concatenate((pads,array),axis=1)

# Rabi frequency calibration
# returns rabi in kHz
@np.vectorize
def volt2rabi(volt):
#     if volt < 0.1 or volt > 5:
#         return 0
    volt = np.log10(volt)
    dbm = 1.5863 +0.2211*volt -0.1022*volt**2 -0.1301*volt**3 -0.0862*volt**4 +0.2323*volt**5 +0.1624*volt**6 -0.1552*volt**7 -0.1206*volt**8
    dbm = 10**dbm
    sqrtpwr = (10**((dbm-30)/10))**(1/2)
    return -0.0332 +0.5832*sqrtpwr -0.0167*sqrtpwr**2

# Provide rabi in kHz
@np.vectorize
def rabi2volt(rabi):
    if rabi <= volt2rabi(0.1) or rabi >= volt2rabi(5):
        print('outside valid range')
        return 0
    def funSolve(v):
        return rabi - volt2rabi(v)
    return brentq(funSolve, 0.1, 5)

# Define an omega32 fit for the tails
def omega32fit(omega,C):
    return (C/np.power(omega,3/2))

# RF spectra class
class RFspectra:
    '''
    A class that analyzes RF spectra data
    '''
    def __init__(self, dataset_name, **kwargs):
        '''Initialize'''
        self.dsname = dataset_name
        # Setup default values for var
        l1 = dict(subsample = 3, cropset = dict(center=(1130, 1600), width=400, height=800), angle=4)
        l2 = dict(Isat=77, time=10, pixel=0.71e-6, sigma=cst.sigma, bg_width=10, od_method='log', bad_light=0)
        l3 = dict(trap_f = 23.9, select = 1, xsec_method = ('poly2',5,0.8), ellipticity=1, center=None, xsec=None)
        l4 = dict(omega0=76.0327, refset=False, nonintset=False, intset=False)
        self.var = {**l1, **l2, **l3,**l4, **kwargs}
        self.omega0 = self.var.get('omega0')

        # Set the directories
        if not os.path.exists('data'):
            os.mkdir('data')
        if not os.path.exists('figures'):
            os.mkdir('figures')


    def set_reference(self, name):
        '''Set a reference image that calibrates the cross-section'''
        self.refname = name
        # Import the reference
        self.hybref = hb.Hybrid(name=self.refname, cropset=self.var.get('cropset'))
        self.hybref.process()
        self.xsec = self.hybref.xsec
        add = dict(refset=True)
        self.var = {**self.var, **add}
        self.zpositions = self.hybref.z

    def inspect_reference(self):
        '''Inspect the reference image'''
        # Check crops
        self.hybref.inspect_crop()
        print('Average Intensity A: {:.3f}'.format(np.nanmean(self.hybref.img.Ii)/(self.hybref.var['Isat']*self.hybref.var['subsample']**2*self.hybref.var['time'])))
        print('Average Count Binned A: {:.3f}'.format(np.nanmean(self.hybref.img.Ii)))

        # Check cross-sections
        fig, ax = plt.subplots(nrows=2, ncols=2,figsize=(8,5))
        self.hybref.xsec.infoplot(axs=[ax[0,0],ax[1,0]],left=self.hybref.denhyb.left_pixel,right=self.hybref.denhyb.right_pixel)

        # Check offset corrections
        im = ax[0,1].imshow(np.transpose(self.hybref.bg.alpha))
        ax[0,1].set_title('Background Gradient')
        cbar_ax = fig.add_axes([0.92, 0.57, 0.02, 0.3])
        fig.colorbar(im, cax=cbar_ax)
        od = np.log(self.hybref.bg.alpha * self.hybref.bg.img.Ii / self.hybref.bg.img.If)
        od_ = np.log(self.hybref.bg.img.Ii / self.hybref.bg.img.If)
        ax[1,1].plot(np.sum(od ,axis=1),'k-')
        ax[1,1].plot(np.sum(od_ ,axis=1),'r--')
        ax[1,1].plot(od*0,'b-')
        ax[1,1].set_title('Offset Correction')

        plt.figure(figsize=(5,3.5))
        plt.plot(self.hybref.z*1e6, self.hybref.n,'.')
        plt.xlabel('z [um]')
        plt.ylabel('n')
        plt.show()

    def set_interacting(self, names):
        '''Analyze a noninteracting spectrum to get the bare resonance'''

        # Set the imagenames
        self.intnames = names
        self.intdata = [hb.Hybrid(name=name, cropset=self.var.get('cropset'), xsec=self.xsec, zpositions=self.zpositions, center=self.hybref.denhyb.center) for name in tqdm(self.intnames, leave=False)]

        # Get the RF freq and power
        tullia = db.Tullia()
        dfout = tullia.image_query(self.intnames, ['RFSpect', 'Spect Volt','PulseTime'])
        self.intrf = dfout['RFSpect'].tolist()[::2]
        self.intrfraw = dfout['RFSpect'].tolist()[::2]
        self.intpower = dfout['Spect Volt'].tolist()[::2]
        self.intrabi = volt2rabi(self.intpower)
        self.pulsetime = 1e-3*dfout['PulseTime'].tolist()[0]

    def interacting_process(self):
        # Process the images
        for hyb in tqdm(self.intdata, leave=False): hyb.process()
        arrivals = self.intdata[0::2]
        refs = self.intdata[1::2]

        # Sort in increasing RF freq
        rfsorter = np.argsort(self.intrf)
        self.rfsorter = rfsorter
        self.intnamessort = [self.intnames[i] for i in rfsorter]
        self.intrf = [self.intrf[i] for i in rfsorter]
        self.intpower = [self.intpower[i] for i in rfsorter]
        refs = [refs[i] for i in rfsorter]
        arrivals = [arrivals[i] for i in rfsorter]

        # Create the main meshgrids and meshgrids for plotting
        self.rfgrid, self.zgrid = np.meshgrid(self.intrf, self.hybref.z)
        self.rfplot = 1e3*(self.rfgrid-self.omega0)  # in kHz
        self.zplot = 1e6*self.zgrid # in um

        # Insert the spectra into 2D arrays
        self.arrivals_spectrum = np.zeros([len(arrivals[0].n),len(arrivals)])
        for i in range(len(arrivals)):
            self.arrivals_spectrum[:,i] = arrivals[i].n
        self.refs_spectrum = np.zeros([len(refs[0].n),len(refs)])
        for i in range(len(refs)):
            self.refs_spectrum[:,i] = refs[i].n

        # Get transfers
        self.transfer = self.arrivals_spectrum/self.refs_spectrum

        # Create a new meshgrid with the detunings normalized to local Fermi energies
        self.EFgrid = 1e-3*cst.n2EFHz(self.refs_spectrum)
        self.kfgrid = np.power(6*(np.pi**2) *self.refs_spectrum,1/3)
        self.nutildeplot = self.rfplot/self.EFgrid

    def get_temperatures(self, tempfit_bins=70, tfit_U_range=[500,3000]):
        self.tempfit_bins = tempfit_bins
        self.tfit_U_range = tfit_U_range
        # Extract the temperatures from the reference shots
        self.thermos = [hbthermo.Thermodynamics(self.zpositions, ref, bins=self.tempfit_bins, tfit_U_range=self.tfit_U_range) for ref in self.refs_spectrum.T]
        self.ttildes = [thermo.Ttilde for thermo in self.thermos]
        self.ttildes = np.array(self.ttildes).T

    def background_subtract(self, bgrefnames, bins=20, spline_smooth=1.3):
        # Subtract backgrounds: no RF arrivals
        self.bgrefnames = bgrefnames
        self.bgrefs = [hb.Hybrid(name=name, cropset=self.var.get('cropset'), xsec=self.xsec) for name in self.bgrefnames]
        for hyb in tqdm(self.bgrefs, leave=False): hyb.process()
        bgs = self.bgrefs[0::2]
        self.bgarray = np.zeros([len(bgs[0].n), len(bgs)])
        for i in range(len(bgs)):
            self.bgarray[:,i] = bgs[i].n
        self.bgprofile = np.nanmean(self.bgarray,1)
        self.bgbin = hbthermo.binxy(1e6*self.zpositions,1e-15*self.bgprofile, bins=bins)
        binx, biny, binerr = self.bgbin
        spl = scipy.interpolate.UnivariateSpline(x=binx,y=biny)
        spl.set_smoothing_factor(spline_smooth)
        self.bgprofile_smooth = 1e15*spl(1e6*self.zpositions)
        self.arrivals_spectrum = self.arrivals_spectrum-np.tile(self.bgprofile_smooth, (np.size(self.arrivals_spectrum,1), 1)).T
        # Get transfers
        self.transfer = self.arrivals_spectrum/self.refs_spectrum

    def background_plot(self):
        # monitor the background subtraction smoothing
        plt.plot(1e6*self.zpositions,self.bgprofile,'b.', label='data')
        plt.plot(1e6*self.zpositions, self.bgprofile_smooth, 'g', lw=3, label='smoothed profile')
        plt.errorbar(self.bgbin[0], 1e15*self.bgbin[1], 1e15*self.bgbin[2], fmt='ro', lw=3, label='binned data' )
        plt.legend()
        plt.show()



    '''Inspection plots'''
    def inspect_arrivals(self):
        figname = 'interacting_arrivals_'+self.dsname+'.png'
        self.plot_array(self.arrivals_spectrum, delta_units='khz', fig_save_name=figname)

    def inspect_refs(self):
        figname = 'interacting_refs_'+self.dsname+'.png'
        self.plot_array(self.refs_spectrum, delta_units='khz', fig_save_name=figname)

    def inspect_transfer(self, clims_input=[0,1], xlims=[0,10]):
        figname = 'interacting_transfer_'+self.dsname+'.png'
        self.plot_array(self.transfer, delta_units='ef', fig_save_name=figname,clims=clims_input, xlims=xlims)

    def inspect_temperatures(self):
        figname = 'temperatures_'+self.dsname+'.png'
        self.plot_array(np.ma.masked_array(self.ttildes,np.logical_or((self.ttildes>2), np.isnan(self.ttildes))), delta_units='khz', fig_save_name=figname)



    '''Contact analysis'''
    def process_contacts(self, method='scaling', OmegaRabi=None, tbin_edges=None, fit_range=None):
        if OmegaRabi==None:
            OmegaRabi=1e3*2*np.pi*self.intrabi[0]
        self.plaincontacts =  self.transfer*self.refs_spectrum*np.power(self.nutildeplot,3/2) *(cst.hbar/cst.mass)* np.power(self.kfgrid,3)
        self.OmegaRabi = OmegaRabi
        self.contact_method = method

        # Set the range of detunings to fit over (if there's data too far in the tails, or outside the linear response regime)
        if fit_range is not None:
            fit_range = 2e3*np.pi*np.array(fit_range)
        self.fit_range = fit_range

        # Process the contacts using the specified method
        if method=='scaling':
            #The scaling method scales the spectra by nutilde^3/2, with some additional factors
            self.contacts = np.sqrt(2)*np.pi*(cst.hbar/cst.mass)*(self.transfer)*np.power(self.nutildeplot,3/2)*(1/(self.pulsetime*self.OmegaRabi**2))*self.kfgrid**2
        elif method=='fitting':
            #The fitting method takes the tails of the spectra and blindly fits omega^3/2
            if tbin_edges==None:
                self.self.tbin_edges = np.array([0.07,0.09,.1,0.11, 0.12,0.13, 0.14, 0.15,0.16, 0.17,0.18,0.19, 0.2, 0.22, .25, .27, 0.3, 0.5, 1, 6])
            else:
                self.tbin_edges = tbin_edges
            self.tbin_edges = np.array(self.tbin_edges)
            # Bin the data by temperatures and fit the tails of each spectrum of constant t/tf
            self.bin_temperatures()
            self.fit_tails()

    '''For contact evaluation from scaled spectra:'''
    def inspect_contacts(self, clims_input=[0,4], xlims=[0,12]):
        figname = 'contacts_'+self.dsname+'.png'
        self.plot_array(self.contacts, delta_units='ef', fig_save_name=figname, clims=clims_input, xlims=xlims)

    def mask_contacts(self, mask_range=[2,4]):
        # mask the contacts
        mask = (self.nutildeplot>mask_range[0])&(self.nutildeplot<mask_range[1])
        self.masked_contacts = ma.masked_array(self.contacts, ~mask)
        self.masked_ttildes = ma.masked_array(self.ttildes, ~mask)

    def inspect_masked_contacts(self, xlims=[0,12]):
        self.plot_array(self.masked_contacts, delta_units='ef', clims=[0,4], xlims=xlims)
        self.plot_array(self.masked_ttildes, delta_units='ef', clims=[0,1], xlims=xlims)

    def bin_contacts(self, bins=140):
        # Compress the masked arrays
        Ttilde, Ctilde = self.masked_ttildes.compressed(), self.masked_contacts.compressed()
        # Mask them
        mask_final = np.isfinite(Ctilde)&np.isfinite(Ttilde)&(Ttilde<3)
        # scatter
        self.Ttilde_scatter, self.Ctilde_scatter = ma.masked_array(Ttilde,~mask_final), ma.masked_array(Ctilde,~mask_final)
        # bin
        list = hbthermo.binxy(self.Ttilde_scatter.compressed(),self.Ctilde_scatter.compressed(),bins=bins)
        self.Ttildebin, self.Ctildebin, self.Ctildeerr = list[0], list[1], list[2]

    def plot_ctilde(self, xlims=[0,1.4], show_virial=True, fig_save_name=''):
        # plot ctilde
        plt.figure(figsize=(5,4))
        plt.errorbar(self.Ttildebin, self.Ctildebin ,yerr=self.Ctildeerr, fmt='o', label='data', capsize=0)
        plt.plot([0.17, 0.17], [0, 10], 'r-', lw=0.5, label='Tc')
        if show_virial:
            vir = virial.VirialUnitarity(BetaMuRange=[-6,0.1])
            CvT = vir.TTilde
            CvC = vir.CI_NkF
            plt.plot(CvT, CvC,'g-', linewidth=2, label='3rd order virial')
        plt.legend()
        plt.xlim(xlims)
        plt.ylim([0,5])
        plt.xlabel(r'$T/T_F$',fontsize=14)
        plt.ylabel(r'$C/Nk_F$',fontsize=14)
        # Save the figure
        if not fig_save_name=='':
            fig_save_name = self.dsname+'_ctilde_scaled'+'.png'
            plt.savefig('figures/'+ fig_save_name, dpi=300)
            print('Figure saved at figures/'+fig_save_name)
        plt.show()


    '''For contact evaluation from fitting the tails of the spectra'''
    def fit_tails(self):
        # Fit the tails and evaluate the contacts
        self.ctilde_fit = np.zeros((self.numttf,1))
        self.ctilde_fit_errs = np.zeros((self.numttf,1))
        self.cfitarray=[]
        self.cfiterrarray=[]
        # For each temperature bin,
        for tposition in range(self.numttf):
            xdata = 2*np.pi*(self.unique_rfs*1e3)
            ydata = self.transfer_over_kf_means[tposition,:]
            intfit = tp.Curve(xdata, ydata)
            fullintfit, err = intfit.fit(omega32fit, (np.nanmean(intfit.y),), xlim = self.fit_range, plot=False);
            self.ctilde_fit[tposition] = 4*np.pi*fullintfit[0]*np.sqrt(cst.mass/cst.hbar)/(self.OmegaRabi**2 *self.pulsetime)
            self.ctilde_fit_errs[tposition] = 4*np.pi*err[0]*np.sqrt(cst.mass/cst.hbar)/(self.OmegaRabi**2 *self.pulsetime)
            self.cfitarray.append(fullintfit)
            self.cfiterrarray.append(err)

        # Flatten the ctilde array, and replace all the zeros (where tbin_edge is outside the temperature range) to nan
        self.ctilde_fit = self.ctilde_fit.flatten()
        self.ctilde_fit[self.ctilde_fit==0]='nan'

    def plot_tail_fits(self, tbin_edges_select=None, fill=True, fig_save_name=''):
        # Plot the omega^-3/2 fits
        if tbin_edges_select==None:
            tbin_edges_select = self.tbin_edges[1::4]

        tbin_edges_select = np.array(tbin_edges_select)
        plt.figure(figsize=(6,6))
        # Select certain T/Tf to plot
        for tbin in tbin_edges_select:
            tposition = np.argmin(np.abs(np.array(self.tbin_edges)-tbin))
            xdata = 2*np.pi*(self.unique_rfs*1e3)
            ydata = self.transfer_over_kf_means[tposition,:]
            errdata = self.transfer_over_kf_errs[tposition,:]
            # Plot the errorbar data
            line = plt.errorbar(self.unique_rfs,ydata,errdata, fmt='o', label=r'$T/T_\mathrm{F} = $'+ '%.2f'%self.tbin_edges[tposition], capsize=0)
            # Plot the fits
            xfit = np.linspace(np.nanmin(xdata), np.nanmax(xdata), 400)
            # Choose whether or not to plot the fit errors as filled solid fit curves
            if fill==False:
                y0 = omega32fit(xfit,*(self.cfitarray[tposition]))
                plt.plot(xfit/(2e3*np.pi),y0, color= line[0].get_color())
            else:
                y1 = omega32fit(xfit,*(self.cfitarray[tposition]+self.cfiterrarray[tposition]))
                y2 = omega32fit(xfit,*(self.cfitarray[tposition]-self.cfiterrarray[tposition]))
                plt.fill_between(xfit/(2e3*np.pi),y1,y2, facecolor=line[0].get_color(),lw=0,alpha=0.5, interpolate=True)

        ax = plt.gca()
        ax.tick_params(labelsize=14)
        plt.legend()
        plt.ylabel(r'$I(\omega)/k_\mathrm{F}$', fontsize=16)
        plt.xlabel(r'$\omega\mathrm{ \;[kHz]}$', fontsize=16)

        # Save the figure
        if not fig_save_name=='':
            fig_save_name = self.dsname+'_tail_fit'+'.png'
            plt.savefig('figures/'+ fig_save_name, dpi=300)
            print('Figure saved at figures/'+fig_save_name)

        plt.show()

    def inspect_binned_arrivals(self):
        self.plot_binned_array(self.arrival_means)
    def inspect_binned_transfers(self):
        self.plot_binned_array(self.transfer_means)
    def inspect_binned_kfs(self):
        self.plot_binned_array(self.kf_means)

    def plot_ctilde_fit(self, xlims=None, show_virial=True, fig_save_name=''):
        # Plot the fitted contacts
        if xlims==None:
            xlims = [0,1]
        plt.figure(figsize=(6,6))
        plt.errorbar(self.tbin_edges, self.ctilde_fit, self.ctilde_fit_errs, fmt='o', label='data', capsize=0)
        plt.plot([0.17, 0.17], [0, 10], 'r-', lw=0.5, label='Tc')
        if show_virial:
            vir = virial.VirialUnitarity(BetaMuRange=[-6,0.1])
            CvT = vir.TTilde
            CvC = vir.CI_NkF
            plt.plot(CvT, CvC,'g-', linewidth=2, label='3rd order virial')
        plt.ylim(0,5)
        plt.xlim(xlims)
        plt.xlabel(r'$T/T_F$',fontsize=14)
        plt.ylabel(r'$C/Nk_F$',fontsize=14)
        plt.legend()

        # Save the figure
        if not fig_save_name=='':
            fig_save_name = self.dsname+'_ctilde_fit'+'.png'
            plt.savefig('figures/'+ fig_save_name, dpi=300)
            print('Figure saved at figures/'+fig_save_name)

        plt.show()


    ''' Code for plotting an array that depends on ttf and detuning '''
    def plot_binned_array(self, array_to_plot, fig_save_name='', clims=None, xlims=None, ylims=None):
        # Plot any array that's a function of rfplot and nutildeplot, no interpolation
        if xlims==None:
            xlims = [np.min(self.rfgrid_bin), np.max(self.rfgrid_bin)]
        if ylims==None:
            ylims = [0.07,0.3]
        xlabel =r'$\delta = (\omega-\omega_0)/2\pi$ [kHz]'

        # Plot the bare arrivals without any warping or normalization
        plt.figure(figsize=(6,6))
        plt.pcolor(self.rfgrid_bin, self.tgrid_bin, array_to_plot)
        plt.xlabel(xlabel, fontsize=14)
        plt.ylabel(r'$T/T_F$', fontsize=14)
        plt.xlim(xlims)
        plt.ylim(ylims)
        plt.colorbar()

        if not clims==None:
            plt.clim(clims)

        # Save the figure
        if not fig_save_name=='':
            plt.savefig('figures/'+fig_save_name+'.png', dpi=300)
            print('Figure saved at figures/'+fig_save_name)

        plt.show()

    '''
    Code for binning the arrivals, transfers, and kfs by T/Tf:
        - run only after temperatures are known
        - currently used for the fitting method of finding contacts
    '''
    def bin_temperatures(self):
        # Bin the data by similar temperatures and equal RF values
        ttildes_positions = np.digitize(self.ttildes, self.tbin_edges)
        self.unique_rfs = np.unique(self.rfplot[1,:])
        self.numrf = len(self.unique_rfs)
        numarrivals = np.shape(self.arrivals_spectrum)[1]
        self.numttf = len(self.tbin_edges)
        self.rfgrid_bin, self.tgrid_bin = np.meshgrid(self.unique_rfs, self.tbin_edges)
        # Initialize the binning arrays
        self.arrival_means = np.zeros((self.numttf,self.numrf))
        self.arrival_errs = np.zeros((self.numttf,self.numrf))
        self.transfer_means = np.zeros((self.numttf,self.numrf))
        self.transfer_errs = np.zeros((self.numttf,self.numrf))
        self.kf_means = np.zeros((self.numttf,self.numrf))
        self.kf_errs = np.zeros((self.numttf,self.numrf))
        self.transfer_over_kf_means = np.zeros((self.numttf,self.numrf))
        self.transfer_over_kf_errs = np.zeros((self.numttf,self.numrf))

        # For each temperature bin,
        for t_bin_index in range(self.numttf):
            # For each RF bin (a unique rf value)
            for rf_bin_index in range(self.numrf):
                arrivals_for_rf_bin = np.array([])
                transfer_for_rf_bin = np.array([])
                kf_for_rf_bin = np.array([])
                transfer_over_kf_for_rf_bin = np.array([])
                # For each arrival profile
                for i in range(numarrivals):
                    # If the rf value of the arrival profile equals the current RF bin
                    if self.rfplot[1,i]==self.unique_rfs[rf_bin_index]:
                        # Find the entries in the current profile that match the temperature bin
                        ttilde_position_column = ttildes_positions[:,1]
                        mask = ttilde_position_column==t_bin_index
                        # Arrivals
                        arrival_column = self.arrivals_spectrum[:,i]
                        masked_arrival_column = np.ma.masked_array(arrival_column,~mask).compressed()
                        arrivals_for_rf_bin = np.concatenate((arrivals_for_rf_bin,masked_arrival_column))
                        # Transfers
                        transfer_column = self.transfer[:,i]
                        masked_transfer_column = np.ma.masked_array(transfer_column,~mask).compressed()
                        transfer_for_rf_bin = np.concatenate((transfer_for_rf_bin,masked_transfer_column))
                        # kfs
                        kf_column = self.kfgrid[:,i]
                        masked_kf_column = np.ma.masked_array(kf_column,~mask).compressed()
                        kf_for_rf_bin = np.concatenate((kf_for_rf_bin,masked_kf_column))
                        # transfer over kfs
                        transfer_over_kf_column = self.transfer[:,i]/self.kfgrid[:,i]
                        masked_transfer_over_kf_column = np.ma.masked_array(transfer_over_kf_column,~mask).compressed()
                        transfer_over_kf_for_rf_bin = np.concatenate((transfer_over_kf_for_rf_bin,masked_transfer_over_kf_column))

                # Throw the arrivals in the appropriate bins
                if len(arrivals_for_rf_bin)!=0:
                    self.arrival_means[t_bin_index, rf_bin_index] = np.nanmean(arrivals_for_rf_bin)
                    self.arrival_errs[t_bin_index, rf_bin_index] = np.nanstd(arrivals_for_rf_bin)/np.sqrt(len(arrivals_for_rf_bin))
                # Throw the transfers in the appropriate bins
                if len(transfer_for_rf_bin)!=0:
                    self.transfer_means[t_bin_index, rf_bin_index] = np.nanmean(transfer_for_rf_bin)
                    self.transfer_errs[t_bin_index, rf_bin_index] = np.nanstd(transfer_for_rf_bin)/np.sqrt(len(transfer_for_rf_bin))
                # Throw the kfs in the appropriate bins
                if len(kf_for_rf_bin)!=0:
                    self.kf_means[t_bin_index, rf_bin_index] = np.nanmean(kf_for_rf_bin)
                    self.kf_errs[t_bin_index, rf_bin_index] = np.nanstd(kf_for_rf_bin)/np.sqrt(len(kf_for_rf_bin))
               # Throw the transfer over kfs in the appropriate bins
                if len(transfer_over_kf_for_rf_bin)!=0:
                    self.transfer_over_kf_means[t_bin_index, rf_bin_index] = np.nanmean(transfer_over_kf_for_rf_bin)
                    self.transfer_over_kf_errs[t_bin_index, rf_bin_index] = np.nanstd(transfer_over_kf_for_rf_bin)/np.sqrt(len(transfer_over_kf_for_rf_bin))




    ''' Code for plotting an array that depends on z and detuning '''
    def plot_array(self, array_to_plot, delta_units='ef', fig_save_name='', clims=None, xlims=None):
        # Plot any array that's a function of rfplot and nutildeplot, no interpolation
        if delta_units=='ef':
            detunings = custom_pad_ef(self.nutildeplot)
            if xlims==None:
                xlims = [-1,8]
            xlabel=r'$\widetilde{\nu}=\delta/E_\mathrm{F}$'
        else:
            detunings = custom_pad_khz(self.rfplot)
            if xlims==None:
                xlims = [np.min(self.rfplot)-5, np.max(self.rfplot)]
            xlabel =r'$\delta = (\omega-\omega_0)/2\pi$ [kHz]'

        # Plot the bare arrivals without any warping or normalization
        plt.figure(figsize=(6,6))
        plt.pcolor(detunings,custom_pad_z(self.zplot),array_to_plot)
        plt.xlabel(xlabel, fontsize=14)
        plt.ylabel(r'z [$\mu$m]', fontsize=14)
        plt.xlim(xlims)
        plt.ylim([np.min(self.zplot), np.max(self.zplot)])
        plt.colorbar()

        if not clims==None:
            plt.clim(clims)

        # Save the figure
        if not fig_save_name=='':
            plt.savefig('figures/'+fig_save_name+'.png', dpi=300)
            print('Figure saved at figures/'+fig_save_name)

        plt.show()

    '''For saving the entire rf object for later use'''
    def save_data(self):
        output = open('data/'+self.dsname+'.pkl', 'wb')
        # Pickle dictionary using protocol 0.
        pickle.dump(self, output, -1)
        output.close()
