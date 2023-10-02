#%%
import numpy as np 
from scipy.special import logit
from astropy.cosmology import Planck18 as cosmo
from numpy.random import default_rng
from wise2mbh.kcorrections import ekcorr, skcorr, lkcorr
rng = default_rng()
""" 
K-correction look-up table from Jarrett+2023 for a set of WISE colors (W1,W1W2,W1W3,W3W4,W2W3)
Input:
morphtype: Morph type of the table to use, only 'E', 'L' and 'S' are supported, if other value is entered, 'E' outputed.

Output:
table: Pandas table of the lookup table for kcorrections, default is 'E' table.
""" 
def kcorr_table(morphtype='E'):
    table = ekcorr
    if morphtype=='S':
        table = skcorr
    elif morphtype=='L':
        table = lkcorr
    return table

"""
Parameter Monte Carlo for error propagation / application of scatter
Inputs:
nominal: Central value of parameter, if scatter, set to 0, default=0
std: Standar deviation of the parameter, if scatter, the scatter itself, default=1
n: size of the random distribution, default=1000

Output:
param: Normal random distribution of size N with mean nominal and standard deviation std
""" 
def param_montecarlo(nominal=0, std=1, n=1000):
    param = nominal + std*rng.standard_normal(number)
    return param

"""
Clipping Dist used to clip W1-W2 color distributions to set values.
Inputs:
array: the W1-W2 color distribution
tresh: limit value to clip
val: default step used to shift, default=0.01
greater_than: If the shift is set to a minimum (True) or a maximum (False), default=True

Output:
new_arrray: W1-W2 color distribution clipped to the tresh value set
""" 
def clipping_dist(array, tresh, val=0.01, greater_than=True):
    mean = np.mean(array, axis=1)
    new_array = np.copy(array)
    for i, m in zip(np.arange(np.shape(array)[0]), mean):
        if greater_than:
            while m > tresh:
                new_array[i,:] -= val
                m = np.mean(new_array[i,:])
        else:
            while m < tresh:
                new_array[i,:] += val
                m = np
    return new_array

"""
Array of Monte Carlo distributions for an array M of quantities
Inputs:
nominal: array-like object of mean values of the quantity
std: error or standard deviation of the values of the quantity
n: size of the Monte Carlo distribution, default=1000

Output:
array_mc: MxN array of random distributions of the quantities. Each row is a normal distribution of mean equal to nominal and standard deviation equal to std
""" 
def array_montecarlo(nominal, std, n=1000):
    to_append = []
    for i,j in zip(nominal,std):
        true_hist = i + j*rng.standard_normal(n)
        to_append.append(true_hist.tolist())

    array_mc = np.array(to_append)
    return array_mc

"""
Distance modulus (mu) for a given redshift
Input:
z: array-like object of redshifts

Output: 
mu: array-like object with the distance modulus
"""
def distance_modulus_z(z):
    if z.size==0:
        return np.array([])  # Return an empty array if z is empty
    lum_dist = cosmo.luminosity_distance(z).value
    mu = 5 * np.log10(lum_dist * 1e6) - 5
    return mu

"""
Luminosity distance for a given redshift
Input:
z: array-like object of redshifts

Output:
lum_dist: array-like object with luminosity distances
"""
def lumdist_z(z):
    if z.size==0:
        return np.array([])  # Return an empty array if z is empty
    lum_dist = cosmo.luminosity_distance(z).value
    return lum_dist

"""
Distance modulus (mu) for a given distance
Input:
dist: array-like object of distances in Mpc

Output:
mu: array-like object with the distance modulus
"""
def distance_modulus_dist(dist):
    if dist.size==0:
        return np.array([])  # Return an empty array if z is empty
    mu = 5*np.log10(dist*1e6)-5
    return mu

"""
WISE Color 2 T-value morphological type.
Inputs:
color: W2-W3 color from WISE bands 

Output:
t_value: Hubble type (T-value) of morphological class of the galaxy
""" 
def w2w3_to_morph(color):
    if type(color) is list:
        color = np.array(color)

    params = [0.745953333333333, 2.7095333333333333, 1.20656567, 1.36157769]
    color_desp_norm = (color-params[0])/params[1]
    t_value = params[2]*logit(color_desp_norm)+params[3]

    t_value = np.where((t_value<-5) | (color_desp_norm<=0), -5, t_value)
    t_value = np.where((t_value>8) | (color_desp_norm>=1), 8, t_value)

    return t_value

"""
WISE to Stellar Mass (M_*) of an extragalactic object
Inputs:
w1abs: array-like object of W1 absolute magnitudes
w1w2: array-like object of W1-W2 colors
resolved: If the WISE mag and colors come from resolved source (like WXSC), it different slopes/intercept, default=False

Output:
log_sm: array-like object with log of the Stellar Mass
""" 
def wise_to_logsm(w1abs, w1w2, resolved=False):
    w1_abs_sun = 3.26                   #from Willmer (2018) https://iopscience.iop.org/article/10.3847/1538-4365/aabfdf

    log_sm_lw1 = -1.96*w1w2 - 0.03
    if resolved:
        log_sm_lw1 = -2.54*w1w2 - 0.17
    log_lw1 = -0.4*(w1abs-w1_abs_sun)   #from Cluver et al. (2014) https://iopscience.iop.org/article/10.1088/0004-637X/782/2/90
    log_sm = log_sm_lw1 + log_lw1
    
    return log_sm

"""
Morphological classification (T_Type) to Bulge-to-Total (B/T)
Input:
t_value: array-like object of morpholofical classification in the Hubble sequence of the extragalactic object

Output:
b_ratio: array-like object of Bulge-to-Total ratios
""" 
def morph_to_bulge_ratio(t_value):
    values = [0.35777546, 7.72843881, 0.09658478, 0.05159308]
    b_ratio =  values[3] + values[0]*(values[1]**(-values[2]*t_value))
    b_ratio[b_ratio>1] = 1

    return b_ratio

"""
Compensated Black Hole Mass (MBH) via aplication of the compensation factor (C_f)
Input:
mbh: array-like object of MBH values

Output:
new_mbh: array-like object of compensated MBH values
""" 
def comp_mbh(mbh):
    params = [-0.18062595,  1.72608017]
    comp_par = params[0]*mbh + params[1]   #derived in notebook for the defined control sample
    new_mbh = mbh + comp_par

    return new_mbh

"""
K-correction of WISE colors
Inputs:
lookup_table: K-correction table to use 
redshift: array-like object of redshifts to calculate K-correction
correction_factor: Color for which to calculate its correction factor, default:W2-W3

Output:
to_return: K-correction factors for selected color
""" 
def get_correction_factor(lookup_table, redshift, correction_factor='W2-W3'):
    if (type(redshift) is float) | (type(redshift) is int):
        if redshift < lookup_table['z'].min() or redshift > lookup_table['z'].max():
            nearest_values = lookup_table.iloc[np.abs(lookup_table['z'] - redshift).argsort()[:2]]
            return np.interp(redshift, nearest_values['z'], nearest_values[correction_factor])
        else:
            return np.interp(redshift, lookup_table['z'], lookup_table[correction_factor])
    else:
        to_return = []
        for redshift in redshift:
            if redshift < lookup_table['z'].min() or redshift > lookup_table['z'].max():
                nearest_values = lookup_table.iloc[np.abs(lookup_table['z'] - redshift).argsort()[:2]]
                a = np.interp(redshift, nearest_values['z'], nearest_values[correction_factor])
            else:
                a = np.interp(redshift, lookup_table['z'], lookup_table[correction_factor])
            to_return.append(a)
        return to_return

"""
K-correction of WISE W1 magnitude
Inputs:
lookup_table: K-correction table to use 
w1: array-like object of W1 magnitudes to calculate K-correction
redshift: array-like object of redshifts to calculate K-correction

Output:
w1_vega_mag_k_corrected: array-like object of K-corrected W1 magnitudes
""" 
def w1_k_corrected(lookup_table,w1,z):
    w1_flux = 309.540 * (10**(-w1/2.5))
    k_fac_w1 = np.array(get_correction_factor(lookup_table=lookup_table, redshift=z, correction_factor='f1'))

    w1_flux_k_corrected = w1_flux*k_fac_w1[:,None]

    w1_ab_mag = -2.5*np.log10(w1_flux_k_corrected) + 8.926
    w1_vega_mag_k_corrected = w1_ab_mag-2.699

    return w1_vega_mag_k_corrected
# %%
