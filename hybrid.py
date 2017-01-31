# General imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import therpy as tp
import scipy
import multiprocessing
from tqdm import tqdm_notebook as tqdm
import pickle
# Specific to interpolation
import os.path
import urllib.request
import numba
import pickle
import math
# For postprocessing
cst = tp.cst(sigmaf=0.5)
from scipy.optimize import curve_fit

# Loading images into memory
class AbsImage:
    def __init__(self, name, bins=(1,1), cropi=None, angle=0, **kwargs):
        # Download Image
        alldata = tp.imageio.imagename2alldata(name)
        # Cropi
        if cropi is None:
            cropi = tp.imagedata.get_cropi(alldata[0], **kwargs)
        # Rotate
        a=angle*math.pi/180.0
        transform=np.array([[np.cos(a),-np.sin(a)],[np.sin(a),np.cos(a)]])
        c_in = np.array(kwargs.get('center'))
        c_out = np.array(kwargs.get('center'))
        offset=c_in-c_out.dot(transform)
        alldata = [scipy.ndimage.interpolation.affine_transform(arr,transform.T,order=2,offset=offset,cval=0.0)for arr in alldata]
        # Subsample and save
        self.If = tp.smooth.subsample2D(alldata[0][cropi], bins)
        self.Ii = tp.smooth.subsample2D(alldata[1][cropi], bins)
        self.name = name

    @property
    def od(self):
        return np.log(self.Ii / self.If)

# Complete fast interpolation code

# Load pre-compiled data
p_ = tp.getpath('Projects','Data','LookupTable','Lookup_Table_Fast_PreCompData_V2.p')
if not os.path.isfile(p_):
    print("Downloading Database -- Might take some time!")
    url = "https://www.dropbox.com/s/4hklnvawtshjay9/Lookup_Table_Fast_PreCompData_V2.p?dl=1"
    u = urllib.request.urlopen(url)
    data = u.read()
    u.close()
    with open(p_, "wb") as f :
        f.write(data)
precompiled_data = pickle.load( open( p_, "rb" ) )

# Jitted interpolation
@numba.jit(nopython=True)
def interp_od_special_jit(IivIn, IfvIn, sim_data):
    # Unload sim_data
    u_si, sf_2d, ocd_2d = sim_data[0], sim_data[1], sim_data[2]
    rows, cols = sf_2d.shape[0], sf_2d.shape[1]
    # Copy inputs and flatten the arrays
    Iiv = IivIn.copy().flatten()  # Flatten so that we can do 1d loop
    Ifv = IfvIn.copy().flatten()  # We will unflatten the arrays when returning
    # Fix low and high OD regions
    bad_low = (Iiv < Ifv)  # For low od (BG), flip Ii and If and make od -> -od
    Iiv[bad_low], Ifv[bad_low] = Ifv[bad_low].copy(), Iiv[bad_low].copy()
    bad_high = (Ifv < 0)   # For high od where If < 0, make If -> -If
    Ifv = np.abs(Ifv)
    # Prepare
    i0v = np.searchsorted(u_si, Iiv)   # Find the indice for closest si
    Pfv = np.zeros_like(Iiv) * np.nan  # Prepare output array, default it with nan
    # Interpolate OD's
    for i in range(Iiv.size):
        Ii, If, i0 = Iiv[i], Ifv[i], i0v[i]
        # Search 4 closest points
        if i0 >= rows or i0 == 0: continue  # If Ii is outside simulation, result is nan
        i1 = np.searchsorted(sf_2d[i0-1,:], If)
        if i1 >= cols: Pfv[i] = 0; continue # If If > max(sf), result is zero atoms
        elif i1 == 0: continue
        i2 = np.searchsorted(sf_2d[i0,:], If)
        if i2 >= cols: Pfv[i] = 0; continue # If If > max(sf), result is zero atoms
        elif i2 == 0: continue
        i0m1 = i0-1
        x1 = u_si[i0m1]
        x2 = u_si[i0]
        dx = x2 - x1
        dx2 = dx**2
        Ary = sf_2d[i0m1, i1-1]
        Bry = sf_2d[i0, i2-1]
        Cry = sf_2d[i0m1, i1]
        Dry = sf_2d[i0, i2]
        Af = ocd_2d[i0m1, i1-1]
        Bf = ocd_2d[i0, i2-1]
        Cf = ocd_2d[i0m1, i1]
        Df = ocd_2d[i0, i2]
        # Interpolate with 4 nearest points
        s = (Ii - x1) / (dx)
        Erx = x1 + (dx) * s
        Ery = Ary + (Bry - Ary) * s
        Frx = x1 + (dx) * s
        Fry = Cry + (Dry - Cry) * s
        Ef = Af + (Bf - Af)  * (((Erx - x1)**2 + (Ery - Ary)**2) / ((dx2 + (Bry - Ary)**2)))**0.5
        Ff = Cf + (Df - Cf)  * (((Frx - x1)**2 + (Fry - Cry)**2) / ((dx2 + (Dry - Cry)**2)))**0.5
        Pfv[i] = Ef + (Ff - Ef) * (((Ii - Erx)**2 + (If - Ery)**2) / (((Frx - Erx)**2 + (Fry - Ery)**2)) )**0.5
    # Make the bad_low od -> -od
    Pfv[bad_low] *= -1
    # Reshape and return
    return Pfv.reshape(*IivIn.shape)

# Wrapper around jitted function to handle passing in pre-compiled data
def interp_od(Ii, If, img_time):
    return interp_od_special_jit(Ii, If, precompiled_data[img_time-1])

# Compute atoms per pixel from raw data

class BackgroundImages:
    def __init__(self,bgod=None,names=None):
        pass

class BorderGradient:
    '''
    2D fit around the edges.
    Computes 2D alpha: alpha*Ii = If
    '''
    def __init__(self, img, width=5):
        self.img = img
        self.width = width
        self.using = self.get_borders()
        self.alpha = self.get_bg_gradient()

    def get_borders(self):
        # Get the border area to be used
        data, w = np.ones_like(self.img.Ii), self.width
        data[w:data.shape[0]-w, w:data.shape[1]-w] = 0
        return data==1

    def get_bg_gradient(self):
        # data = If/Ii
        data = self.img.If / self.img.Ii
        # Setup grid
        xx, yy = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
        # Fit linear gradient
        xy_ = (xx[self.using], yy[self.using])
        z_ = data[self.using]
        good = np.isfinite(z_)
        xy_ = (xy_[0][good], xy_[1][good])
        z_ = z_[good]
        guess = (1, 1e-3, 1e-3)
        fitres, fiterr = curve_fit(self.polynomial2D, xy_, z_, p0=guess)
        return self.polynomial2D_separate(xx, yy, *fitres)

    def infoplots(self):
        plt.figure()
        plt.imshow(self.alpha)
        plt.colorbar()
        plt.figure()
        od = np.log(self.alpha * self.img.Ii / self.img.If)
        od_ = np.log(self.img.Ii / self.img.If)
        plt.plot(np.sum(od ,axis=1),'k-')
        plt.plot(np.sum(od_ ,axis=1),'r--')
        plt.plot(od*0,'b-')

    def polynomial2D(self, xy, b, m1, m2):
        return b + m1*xy[0] + m2*xy[1]

    def polynomial2D_separate(self, x, y, b, m1, m2):
        return b + m1*x + m2*y

class AtomNum:
    def __init__(self, img, Nsat, bg_width, time, pixel, sigma, method='raw', bad_light=0):
        self.bg = BorderGradient(img, bg_width)
        self.od_ = np.log(self.bg.alpha * img.Ii/ img.If) + (self.bg.alpha * img.Ii - img.If)/Nsat
        # process inputs
        if method == 'table':
            self.si = self.bg.alpha * img.Ii/Nsat * (1-bad_light)
            self.sf = img.If/Nsat - bad_light * self.si
            self.od = interp_od(self.si, self.sf, time)
        elif method == 'raw':
            #print('raw, i.e., using log + dI')
            self.si = self.bg.alpha * img.Ii/Nsat * (1-bad_light)
            self.sf = img.If/Nsat - bad_light * self.si
            self.od = np.log(self.si / self.sf) + (self.si - self.sf)
        elif method == 'log':
            #print('log only')
            self.si = self.bg.alpha * img.Ii/Nsat * (1-bad_light)
            self.sf = img.If/Nsat - bad_light * self.si
            self.od = np.log(self.si / self.sf)
        self.atoms = self.od / sigma * pixel**2




@np.vectorize
def area_partial_ellipse(A, a, b=None):
    if b is None: b = a
    if A >= a: return np.pi*a*b
    return 2*A*b*np.sqrt(1-(A/a)**2) + 2*a*b*np.arcsin(A/a)

class XSectionHybrid:
    '''
    Inputs:
        1) data -- OD like
        2) ellipticity
        3) method = ('linear', slice_wdith, fit_range)

    Procedure:
        1) Approximate Center
        2) Fit circles until height drops to xx%, store left, right, center
        3) Extrapolate
        4) Functions for left, right, center, radius, area, sub_area(l, r)
    '''
    def __init__(self, data, ellipticity=1, method=[None]*3):
        # Process Inputs
        self.data = data
        self.ellipticity = ellipticity
        self.method = method
        # Get fitting range
        self.z_edges, self.z_center = self.circle_fitting_range()
        # Fit circles
        self.center_fit, self.radius_fit = self.fit_circles()
        # Extrapolate
        self.z, self.center, self.radius = self.extrapolate()

    def circle_fitting_range(self):
        '''
        Get approximate center and radius by a Thomas-Fermi fit
        Use region radius*fit_range for circle fitting
        '''
        # Inputs
        slice_width = 5 if self.method[1] is None else self.method[1]
        stop_amplitude = 0.75 if self.method[2] is None else self.method[2]
        # Integrate hybrid
        c = tp.Curve(y = np.nanmean(self.data, axis=1))
        c.removenan()
        # Fit HarmonicTF
        def fitfun(x, x0, rad, amp, m, b):
            y = np.real(amp * (1 - ((x-x0)/(rad))**2)**(3/2))
            y[np.isnan(y)] = 0
            y += m*x + b
            return y
        guess = (c.x.shape[0]/2, c.x.shape[0]/4, np.max(c.y), np.max(c.y)/100, np.max(c.y)/1000)
        fitres = c.fit(fitfun, guess, plot=False)[0]
        center, radius_use = round(fitres[0]), round(fitres[1] * stop_amplitude)
        # z_edges and z_center arrays
        z_edges = np.arange(center - radius_use, center + radius_use, slice_width, dtype=np.int)
        z_center = z_edges[0:-1] + (z_edges[1]-z_edges[0])/2.0 - 0.5
        return (z_edges, z_center)

    def fit_circles(self):
        '''
        Fit circles to the range specified
        Measure center and radius at each point
        '''
        # Inputs
        z_center, z_edges = self.z_center, self.z_edges
        # Prepare arrays
        center, radius = np.zeros_like(z_center), np.zeros_like(z_center)
        # Fit circles to each slices
        use_data = self.data.copy()
        use_data[~np.isfinite(use_data)] = 0
        for i in range(self.z_center.size):
            c = tp.Curve(y = np.nanmean(use_data[z_edges[i]:z_edges[i+1],:], axis=0))
            c.removenan()
            guess = (c.x.shape[0] / 2, c.x.shape[0] / 4, np.max(c.y), np.max(c.y)/10, np.max(c.y)/100)
            fitres = c.fit(self.fitfun_circle, guess, plot=False)[0]
            center[i], radius[i] = fitres[0], fitres[1]
        # return results
        return (center, radius)

    def extrapolate(self):
        '''
        Extrapolate the fitted center and radius
        using either polyN or splineN method
        '''
        # Inputs
        method = 'poly4' if self.method[0] is None else self.method[0]
        z_center_fit, center_fit, radius_fit = self.z_center, self.center_fit, self.radius_fit
        # Empty arrays for storage
        z, center, radius = np.arange(self.data.shape[0]), np.arange(self.data.shape[0]), np.arange(self.data.shape[0])
        # Linearly extend the center
        fitres = np.poly1d(np.polyfit(z_center_fit, center_fit, deg=1))
        center = fitres(z)
        # polyN
        if method[0:4] == 'poly':
            fitres = np.poly1d(np.polyfit(z_center_fit, radius_fit, deg=int(method[4:])))
            radius = fitres(z)
            radius[z<z_center_fit[0]] = fitres(z_center_fit[0])
            radius[z>z_center_fit[-1]] = fitres(z_center_fit[-1])
        # splineN
        elif method[0:6] == 'spline':
            fitres = scipy.interpolate.splrep(z_center_fit, radius_fit, s=int(method[6:]))
            radius = scipy.interpolate.splev(z, fitres, der=0)
            radius[z<z_center_fit[0]] = scipy.interpolate.splev(z_center_fit[0], fitres, der=0)
            radius[z>z_center_fit[-1]] = scipy.interpolate.splev(z_center_fit[-1], fitres, der=0)
        # Return
        return (z, center, radius)

    '''
    Useful calls to get center, radius, left, right, area, and sub_area
    '''
    def get_center(self, z):
        z = np.array(z, dtype=np.int32)
        return self.center[z]

    def get_radius(self, z):
        z = np.array(z, dtype=np.int32)
        return self.radius[z]

    def get_left(self, z):
        return self.get_center(z) - self.get_radius(z)

    def get_right(self, z):
        return self.get_center(z) + self.get_radius(z)

    def get_area(self, z):
        return np.pi * self.get_radius(z)**2 * self.ellipticity

    def get_subarea(self, z, l, r):
        a = self.get_radius(z)
        b = a * self.ellipticity
        Al = self.get_center(z) - l
        Ar = r - self.get_center(z)

        # Check for errors
        if np.any(Al <= 0) or np.any(Ar <= 0):
            print("Illegar left and right points given to XSectionHybrid.get_subarea. Returned total area.")
            return self.get_area(z)

        return area_partial_ellipse(Al,a,b)/2 + area_partial_ellipse(Ar,a,b)/2


    def infoplot(self, axs=None, left=None, right=None):
        '''
        Useful information plots: data with fitted center and radius + extrapolation
        Ability to plot on provided axes
        '''
        if axs is None:
            fig, axs = plt.subplots(figsize=(5,5), nrows=2)
        axs[0].imshow(self.data.T, cmap='viridis', aspect='auto', origin='lower')
        axs[0].plot(self.z, self.center,'w--',alpha=0.5)
        axs[0].plot(self.z, self.center - self.radius,'w--',alpha=0.5)
        axs[0].plot(self.z, self.center + self.radius,'w--',alpha=0.5)
        axs[0].scatter(self.z_center, self.center_fit - self.radius_fit,color='white', s=2)
        axs[0].scatter(self.z_center, self.center_fit + self.radius_fit,color='white', s=2)
        axs[0].scatter(self.z_center, self.center_fit,color='white', s=2)
        axs[0].set(xlim=(self.z[0],self.z[-1]))
        axs[0].set_axis_off()

        if left is not None and right is not None:
            axs[0].plot(left,'r-',alpha=0.7)
            axs[0].plot(right,'r-',alpha=0.7)

        axs[1].scatter(self.z_center, self.radius_fit,color='red')
        axs[1].plot(self.z, self.radius,'k')
        axs[1].set(xlim=(self.z[0],self.z[-1]))

    def fitfun_circle(self, x, x0, rad, amp, m, b):
        y = 1 - ((x - x0) / rad) ** 2
        y[y <= 0] = 0
        y[y > 0] = np.sqrt(y[y > 0]) * amp
        y += m * x + b
        return y


class DensityHybrid:
    '''
    Converts 2D atoms/pixel to 1D n(z).
    Allows for integrating entire hybrid, or select partial region.
    '''
    def __init__(self, atoms, xsec, pixel, trap_f=23.9, select=1, center=None):
        # Process Inputs
        self.atoms = atoms
        self.xsec = xsec
        self.pixel = pixel
        self.trap_w = 2*np.pi*trap_f
        self.select = select
        self.center = center
        # Compute n(pixel)
        self.i, self.n = self.density_pixel()
        # Compute n(z) and n(u)
        self.z, self.u = self.get_position()

    def density_pixel(self):
        '''
        integrate over all region, or select partial
        '''
        n = np.nansum(self.atoms, axis=1)
        i = np.arange(n.size)
        self.left_pixel = i*0
        self.right_pixel = i*0 + self.atoms.shape[1] - 1
        # If select = 1 or 0, return entire integration
        if self.select == 1 or self.select == 0:
            n /= self.xsec.get_area(i) * self.pixel**3
            return (i, n)
        # Extract center and radius
        center = self.xsec.get_center(i)
        radius = self.xsec.get_radius(i)
        # If select < 1, use +- round(radius * select)
        if self.select < 1:
            l = np.array(np.round(center-radius*self.select),dtype=np.int)
            r = np.array(np.round(center+radius*self.select),dtype=np.int)
            for j in range(n.size):
                if l[j]<(center[j]-radius[j]): l[j] = round(center[j]-radius[j])
                if r[j]>(center[j]+radius[j]): r[j] = round(center[j]+radius[j])
                n[j] = np.nansum(self.atoms[j,l[j]:r[j]+1])
            n /= self.xsec.get_subarea(i, l, r) * self.pixel**3
            self.left_pixel = l
            self.right_pixel = r
            return (i, n)
        # If select > 1, use fixed +- round(select) pixels
        if self.select > 1:
            l = np.array(np.round(center-self.select),dtype=np.int)
            r = np.array(np.round(center+self.select),dtype=np.int)
            for j in range(n.size):
                if l[j]<(center[j]-radius[j]): l[j] = round(center[j]-radius[j])
                if r[j]>(center[j]+radius[j]): r[j] = round(center[j]+radius[j])
                n[j] = np.nansum(self.atoms[j,l[j]:r[j]+1])
            n /= self.xsec.get_subarea(i, l, r) * self.pixel**3
            self.left_pixel = l
            self.right_pixel = r
            return (i, n)

    def get_position(self):
        '''
        find center and convert pixel to position and potential
        '''
        # Find trap center if not provided
        if self.center is None:
            c = tp.Curve(self.i, self.n)
            c.removenan()
            c = c.copy(y = c.y / np.max(c.y))
            #guess=(c.x.size/2, 1, 1e-8, 1e-4, 1e-8, 1e-12, 1e-16, 1e-20)
            #center = c.fit(self.fitfun_evenpoly, guess, plot=False)[0][0]
            guess = (c.x.size/2, c.x.size/4, c.maxy, c.maxy/10, c.maxy/100)
            self.center = c.fit(self.fitfun_TF, guess, plot=False)[0][0]
        # Find z and u
        z = (self.i - self.center) * self.pixel
        u = 0.5 * cst.mass * self.trap_w**2 * z**2
        return (z, u)

    def fitfun_evenpoly(self, x, x0, b, m, a2=0, a4=0, a6=0, a8=0, a10=0):
        return b + m*(x-x0) + a2*(x-x0)**2 + a4*(x-x0)**4 + a6*(x-x0)**6 + a8*(x-x0)**8 + a10*(x-x0)**10

    def fitfun_TF(self, x, x0, rad, amp, m, b):
        y = np.real(amp * (1 - ((x-x0)/(rad))**2)**(3/2))
        y[np.isnan(y)] = 0
        y += m*x + b
        return y


class Hybrid:
    '''
    A Class that creates the data object from an image name, given known parameters
    '''
    def __init__(self, name, **kwargs):
        # Initialize
        self.name = name
        # Setup default values for var
        l1 = dict(subsample = 3, cropset = dict(center=(1130, 1600), width=400, height=800), angle=4)
        l2 = dict(Isat=77, time=10, pixel=0.6e-6, sigma=cst.sigma, bg_width=10, od_method='table', bad_light=0)
        l3 = dict(trap_f = 23.9, select = 1, xsec_method = ('poly2',10,0.9), ellipticity=1, center=None, xsec=None)
        self.var = {**l1, **l2, **l3, **kwargs}

    def inspect_crop(self):
        # Inspect the crop
        tp.AbsImage(self.name).cropi(**self.var.get('cropset'), plot=True)

    def inspect_atoms(self):
        plt.imshow(self.atoms)
        plt.colorbar()
        plt.show()

    def process(self,**kwargs):
        # update self.var
        self.var = {**self.var, **kwargs}
        self.load_image(subsample=self.var.get('subsample'),
                        angle=self.var.get('angle'),
                        **self.var.get('cropset'))
        self.compute_od(Isat = self.var.get('Isat'),
                        time = self.var.get('time'),
                        pixel= self.var.get('pixel'),
                        sigma= self.var.get('sigma'),
                        bg_width = self.var.get('bg_width'),
                        od_method = self.var.get('od_method'),
                        bad_light = self.var.get('bad_light'))
        self.compute_density(trap_f=self.var.get('trap_f'),
                             select = self.var.get('select'),
                             xsec_method = self.var.get('xsec_method'),
                             pixel = self.var.get('pixel'),
                             ellipticity = self.var.get('ellipticity'),
                             center = self.var.get('center'),
                             xsec = self.var.get('xsec'))

    def load_image(self, subsample=1, angle=0, **cropset):
        '''
        Inputs:
            subsample : determines the subsampling of the image
            **cropset : dict containing one or more of center, width, height, point1, point2, cropi
        Procedure:
            Download Image ==> Crop ==> Subsample ==> Store
        Results:
            self.img  ==  Object with Ii and If properties
            self.var  ==  update with bins, cropset
        '''
        # Download Image
        self.img = AbsImage(name=self.name, bins=(subsample,subsample), angle=angle, **cropset)
        # Update self.var
        self.var['subsample'] = subsample
        self.var['cropset'] = cropset

    def compute_od(self, Isat, time, pixel, sigma, bg_width, od_method, bad_light):
        '''
        Inputs:
            Isat      : saturation counts/us. It will be scaled autometically to account for subsampling and img_time
            time      : imaging time used for LookupTable
            pixel     : pixel size onto atoms. It will be scaled autometically to account for subsampling
            sigma     : effective cross section
            bg_width  : width for border gradient background fix
            od_method : [str]
            bad_light : fraction of Ii that is wrong polarization. Usually 0.03 ish
        Procedure:

            Fix background gradients ==> Compute Ii/Isat, If/Isat ==> Compute OD ==> Fix nan
        Results:
            self.od      : true optical density
            self.atoms   : atoms per pixel
            self.atomnum : AtomNum Object
            self.var     : Isat, time, pixel, sigma, bg_width, od_method, bad_light
        '''
        # Atom Number calculations
        atomnum = AtomNum(img = self.img, Nsat = Isat*self.var['subsample']**2*time,
                               bg_width = bg_width, time = time,
                               pixel = pixel*self.var['subsample'], sigma = sigma,
                               method = od_method,
                               bad_light = bad_light)
        self.od = atomnum.od
        self.atoms = atomnum.atoms
        self.bg = atomnum.bg
        # Update self.var
        add = dict(Isat=Isat, time=time, pixel=pixel, sigma=sigma, bg_width=bg_width, od_method=od_method, bad_light=bad_light)
        self.var = {**self.var, **add}

    def compute_density(self, trap_f, select, xsec_method, pixel, ellipticity, center, xsec):
        '''
        Inputs:
            trap_f       :  trapping frequency in Hz (or provide trap_w)
            xsec_method  :  (str extension, int slice_width, float fit_range)
            select       :  region to use for density calculation: 0,1 for all or number of pixels 2,3,...
            ellipticity  :  1
            center       :  None, center of trap in pixel to force
        Procedure:
            Get XSection ==> Density in Selected Region ==> Center ==> n(z) and n(u)
        Results:
            self.xsec    : XSectionHybrid Object
            self.denhyb  : DensityHybrid Object
            self.n, z, u : arrays for n, z, u of equal length
            self.var     : trap_f, select, xsec_method
        '''
        # Density calculations
        self.xsec = self.var.get('xsec')
        if self.xsec is None:
            self.xsec = XSectionHybrid(data = self.od, method = xsec_method, ellipticity=ellipticity)
        self.denhyb = DensityHybrid(atoms= self.atoms, xsec = self.xsec, select = select,
                                    trap_f = trap_f, pixel=pixel*self.var['subsample'],
                                    center = center)
        self.n, self.z, self.u = self.denhyb.n, self.denhyb.z, self.denhyb.u
        # Update self.var
        add = dict(trap_f=trap_f, select=select, xsec_method=xsec_method, center=center, ellipticity=ellipticity, xsec=self.xsec)
        self.var = {**self.var, **add}
