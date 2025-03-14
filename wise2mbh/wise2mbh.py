import numpy as np 
from scipy.special import logit
from astropy.cosmology import Planck18 as cosmo
from numpy.random import default_rng
from wise2mbh.kcorrections import ekcorr, skcorr, lkcorr
rng = default_rng()

mc_size = 1000

def kcorr_table(morphtype='E'):
    """ 
    K-correction look-up table from Jarrett+2023 for a set of WISE colors (W1,W1W2,W1W3,W3W4,W2W3).

    Input:
    morphtype (str, optional): Morph type of the table to use, only 'E', 'L' and 'S' are supported.

    Output:
    table (pandas.DataFrame): Lookup table for kcorrections, default is 'E' table.
    """ 
    table = ekcorr
    if morphtype=='S':
        table = skcorr
    elif morphtype=='L':
        table = lkcorr
    return table

def param_montecarlo(nominal=0, std=1, n=mc_size):
    """
    Parameter Monte Carlo for error propagation / application of scatter.

    Inputs:
    nominal (float, optional): Central value of parameter, if scatter, set to 0
    std (float, optinal): Standar deviation of the parameter, if scatter, the scatter itself
    n (int, optinal): size of the random distribution

    Output:
    param (numpy.ndarray): Normal random distribution of size N with mean nominal and standard deviation std
    """ 
    param = nominal + std*rng.standard_normal(n)
    return param

def clipping_dist(array, tresh, val=0.01, greater_than=True):
    """
    Clipping Dist used to clip W1-W2 color distributions to set values.

    Inputs:
    array (numpy.ndarray): the W1-W2 color distribution
    tresh (float): limit value to clip
    val (float, optional): default step used to shift
    greater_than (boolean, optional): If the shift is set to a minimum (True) or a maximum (False)

    Output:
    new_arrray (numpy.ndarray): W1-W2 color distribution clipped to the tresh value set
    """ 
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
                m = np.mean(new_array[i,:])
    return new_array

def array_montecarlo(nominal, std, n=mc_size):
    """
    Array of Monte Carlo distributions for an array M of quantities

    Inputs:
    nominal (numpy.ndarray): array-like object of mean values of the quantity
    std (numpy.ndarray): error or standard deviation of the values of the quantity
    n (int, optional): size of the Monte Carlo distribution

    Output:
    array_mc (numpy.ndarray): MxN array of random distributions of the quantities. Each row is a normal distribution of mean equal to nominal and standard deviation equal to std
    """ 
    to_append = []
    for i,j in zip(nominal,std):
        true_hist = i + j*rng.standard_normal(n)
        to_append.append(true_hist.tolist())

    array_mc = np.array(to_append)
    return array_mc

def distance_modulus_z(z):
    """
    Distance modulus (mu) for a given redshift

    Input:
    z (numpy.ndarray): array-like object of redshifts

    Output: 
    mu (numpy.ndarray): array-like object with the distance modulus
    """
    if z.size==0:
        return np.array([])  # Return an empty array if z is empty
    lum_dist = cosmo.luminosity_distance(z).value
    mu = 5 * np.log10(lum_dist * 1e6) - 5
    return mu

def lumdist_z(z):
    """
    Luminosity distance for a given redshift

    Input:
    z (numpy.ndarray): array-like object of redshifts

    Output:
    lum_dist (numpy.ndarray): array-like object with luminosity distances
    """
    if z.size==0:
        return np.array([])  # Return an empty array if z is empty
    lum_dist = cosmo.luminosity_distance(z).value
    return lum_dist

def distance_modulus_dist(dist):
    """
    Distance modulus (mu) for a given distance

    Input:
    dist (numpy.ndarray): array-like object of distances in Mpc

    Output:
    mu (numpy.ndarray): array-like object with the distance modulus
    """
    if dist.size==0:
        return np.array([])  # Return an empty array if z is empty
    mu = 5*np.log10(dist*1e6)-5
    return mu

def w2w3_to_morph(color, mc=False, n=mc_size):
    """
    WISE Color 2 T-value morphological type

    Inputs:
    color (numpy.ndarray): W2-W3 color from WISE bands 
    mc (boolean, optinal): If the calculations consider Monte Carlo approach

    Output:
    t_value: Hubble type (T-value) of morphological class of the galaxy
    """ 
    if type(color) is list:
        color = np.array(color)

    params = [0.745953333333333, 2.7095333333333333, 1.20656567, 1.36157769]
    color_desp_norm = (color-params[0])/params[1]
    t_value = params[2]*logit(color_desp_norm)+params[3]
    if mc:
        params_mc = [1.20656567, 1.36157769]
        t_value = param_montecarlo(params_mc[0],0.01,n=n)*logit(color_desp_norm)+param_montecarlo(params_mc[1],0.02,n=n)

    t_value = np.where((t_value<-5) | (color_desp_norm<=0), -5, t_value)
    t_value = np.where((t_value>8) | (color_desp_norm>=1), 8, t_value)

    return t_value

def wise_to_logsm(w1abs, w1w2, resolved=False):
    """
    WISE to Stellar Mass (M_*) of an extragalactic object

    Inputs:
    w1abs (numpy.ndarray): array-like object of W1 absolute magnitudes
    w1w2 (numpy.ndarray): array-like object of W1-W2 colors
    resolved (boolean, optional): If the WISE mag and colors come from resolved source (like WXSC), it uses different slopes/intercept

    Output:
    log_sm (numpy.ndarray): array-like object with log of the Stellar Mass
    """ 
    w1_abs_sun = 3.26                   #from Willmer (2018) https://iopscience.iop.org/article/10.3847/1538-4365/aabfdf

    log_sm_lw1 = -1.96*w1w2 - 0.03
    if resolved:
        log_sm_lw1 = -2.54*w1w2 - 0.17
    log_lw1 = -0.4*(w1abs-w1_abs_sun)   #from Cluver et al. (2014) https://iopscience.iop.org/article/10.1088/0004-637X/782/2/90
    log_sm = log_sm_lw1 + log_lw1
    
    return log_sm

def morph_to_bulge_ratio(t_value):
    """
    Morphological classification (T_Type) to Bulge-to-Total (B/T)

    Input:
    t_value (numpy.ndarray): array-like object of morpholofical classification in the Hubble sequence of the extragalactic object

    Output:
    b_ratio (numpy.ndarray): array-like object of Bulge-to-Total ratios
    """ 
    values = [0.35777546, 7.72843881, 0.09658478, 0.05159308]
    b_ratio =  values[3] + values[0]*(values[1]**(-values[2]*t_value))
    b_ratio[b_ratio>1] = 1

    return b_ratio

def bulge_to_mbh(bulgemass, mc=False, n=mc_size):
    """
    Bulge Mass to Black Hole Mass from Schutte+2019

    Input:
    bulgemass (numpy.ndarray): array-like object of bulge masses
    mc (boolean, optional): Boolean to use monte carlo parameters, default=False

    Output:
    mbh (numpy.ndarray): array-like object of Black hole masses
    """ 
    mbh = 1.24*(bulgemass - 11) +8.8

    if mc:
        param1 = param_montecarlo(1.24,0.08, n=n)
        param2 = param_montecarlo(8.8,0.09, n=n)
        sct = param_montecarlo(0,0.68, n=n)

        mbh = param1*(bulgemass - 11) + param2 + sct
    return mbh

def comp_mbh(mbh):
    """
    Compensated Black Hole Mass (MBH) via aplication of the compensation factor (C_f)

    Input:
    mbh (numpy.ndarray): array-like object of MBH values

    Output:
    new_mbh (numpy.ndarray): array-like object of compensated MBH values
    """ 
    params = [-0.10378699,  0.98153064]
    comp_par = params[0]*mbh + params[1]   #derived in notebook for the defined control sample
    new_mbh = mbh + comp_par

    return new_mbh

def get_correction_factor(lookup_table, redshift, correction_factor='W2-W3'):
    """
    K-correction of WISE colors

    Inputs:
    lookup_table (pandas.DataFrame): K-correction table to use 
    redshift (float, int or numpy.ndarray): array-like object of redshifts to calculate K-correction
    correction_factor (str, optional): Color for which to calculate its correction factor, default: 'W2-W3'

    Output:
    to_return (float or numpy.ndarray): K-correction factors for selected color
    """ 
    if isinstance(redshift, (float,int,np.float64,np.float32)):
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

def w1w2_treshold_qso(w2w3):
    """
    W1-W2 limit in boxy region

    Inputs:
    w2w3 (numpy.ndarray): Observed W2-W3 color for a set of sources to test

    Output:
    w1w2_tresh (numpy.ndarray): Boxy limit in W1-W2 color
    """ 
    w1w2_tresh = (0.05*w2w3)+0.38
    return w1w2_tresh

def w3_to_SFR(w3,w2w3,z,ulirgs=False,mc=False,n=mc_size):
    """
    log SFR from W3 magnitude, using Cluver+17 relation

    Inputs:
    w3 (float or numpy.ndarray): W3 magnitude, in mag
    w2w3 (float or numpy.ndarray): Observed W2-W3 color
    z (float): Redshift of the object
    mc (bool): Use of MC error propagation, default: False
    n (int): Size of the MC array used, default: mc_size

    Output:
    wise_sfr (float or numpy.ndarray): log SFR (in M_Sun/yr)
    """ 
    f_v0 = 29.045                       #zero magnitude flux density fora assumed constant power-law (https://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html)
    sun_ergs = np.log10(3.839e33)       #from Jarret+13 in ergs/s (https://iopscience.iop.org/article/10.1088/0004-6256/145/1/6)
    flux_dens = f_v0*(10**(- w3/2.5))   #Eq. 1 from same documentation as above, in Jansky (10-26 W m-2 Hz-1)

    consts = [0.02172314, 2.09222101, 0.91801391]
    fc = consts[0] * ((w2w3-consts[1])**2) + consts[2]  #get correction factor based on W2-W3 obs color from Wright+10 (https://iopscience.iop.org/article/10.1088/0004-6256/140/6/1868/pdf)
    flux_dens_corr = flux_dens/fc                       #flux corrected by color

    hz_W3 = 3e8/1.2e-5                  #12 micrometers to Hz (for vLv)

    dist_in_mpc = cosmo.luminosity_distance(z).value
    dist_in_m = dist_in_mpc*3.085e22    # to meters

    flux_watts = flux_dens_corr*10e-26          # in W m-2 Hz-1
    lum_nu = 4*np.pi*flux_watts*(dist_in_m**2)  # in W Hz-1
    nu_lum_nu = hz_W3*lum_nu                    # in W

    nu_lum_nu_ergs = nu_lum_nu*1e6              # W to ergs/s
    log_nu_lum = np.log10(nu_lum_nu_ergs)       # to log scales
    log_nu_lum_lsun = log_nu_lum - sun_ergs     # divided by solar luminosity

    wise_sfr = 0.889*log_nu_lum_lsun - 7.76     # from Cluver+17 (https://iopscience.iop.org/article/10.3847/1538-4357/aa92c7/pdf)

    if ulirgs:
        wise_sfr = 1.43*log_nu_lum_lsun - 13.17 # from Cluver+17 but for AGN/ULIRGS (https://iopscience.iop.org/article/10.3847/1538-4357/aa92c7/pdf)

    if mc & (not ulirgs):
        param1 = param_montecarlo(0.889,0.018,n=n)
        param2 = param_montecarlo(-7.76,0.15,n=n)
        scatter = param_montecarlo(0,0.15,n=n)

        wise_sfr = param1*log_nu_lum_lsun + param2 + scatter
    
    if mc & ulirgs:
        param1 = param_montecarlo(0.889,0.018,n=n)
        param2 = param_montecarlo(-7.76,0.15,n=n)
        scatter = param_montecarlo(0,0.15,n=n)

        wise_sfr = param1*log_nu_lum_lsun + param2 + scatter

        ids = np.where(log_nu_lum_lsun.mean(axis=1)>=10)[0]

        param1 = param_montecarlo(1.43,0.161,n=n)
        param2 = param_montecarlo(-13.17,1.66,n=n)
        scatter = param_montecarlo(0,0.22,n=n)

        wise_sfr_ulirgs = param1*log_nu_lum_lsun[ids] + param2 + scatter

        wise_sfr[ids] = wise_sfr_ulirgs

    return wise_sfr

def agn_fraction(w1w2):
    """
    AGN fraction to consider to avoid contamination from Hernandez-Yevenes+25

    Inputs:
    w1w2 (float or numpy.ndarray): Observed W1-W2 color

    Output:
    agn_frac (float or numpy.ndarray): Fraction of the total Mstellar to consider for AGN contaminated sources
    """ 
    agn_frac = 0.12 - 0.26 * np.log10(w1w2)
    return agn_frac

# FUNCTIONS FOR WXSC TUTORIALS ONLY, NOT IMPORTANT AND VERY SPECIFIC

def drop_irregulars(strings):
    def criteria(string):
        for i in range(len(string)-1):
            if ((string[i].isupper()) and (string[i + 1] == 'm')) | (string[0]=='I'):
                return False
        return True
    boolean_array = [criteria(string) for string in strings]
    return np.array(boolean_array)

def keep_first_and_non_capitals(strings):
    def criteria(string):
        result = string[0]  # Initialize the result with the first letter

        # Iterate through the remaining characters in the string
        for char in string[1:]:
            if char.islower() or not char.isalpha():
                result += char

        return result.strip()
    output = [criteria(s) for s in strings]
    return np.array(output)