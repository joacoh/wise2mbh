'''
WISE2MBH Internal Pipeline

This is a modified version of the WISE2MBH internal pipeline for the ETHER database. 
This pipeline may also be outdated compared to the public version, and it will not be updated.
Only non-vulnerable processes are shown in this pipeline and this version does not work in other systems. 
If you want to use a functional pipeline, please use the 'pipeline.py' file shared in the repo.

ONLY FOR TRANSPARENCY, THIS PIPELINE DOES NOT WORK.
'''
import wise2mbh as wm 
import numpy as np
from astropy.table import Table, unique
import time

print('Welcome to WISE2MBH Internal, now uploading and depuring samples...')

no_data = 9876543.0
chunk = 2000
mc_size = int(1e3)

print('Running with an MC size of {}, if this is not correct, please abort.'.format(mc_size))

'''
Importing and Depuring initial sample
'''

input_sample = Table.read(directory)
ether = Table.read(ether_directory)

input_sample.sort(keys='angDist')

print('Cleaning old estimates from ETHER subsample ({})...'.format(ether_directory))

to_clean_id = [idx for idx, val in enumerate(ether['MBHPUB']) if val == 50]

columns_to_clean_1 = ['MBH','MBHLO','MBHHI','MBHSOURCE','MBHPUB','MBHLIM','MBHQUALITY']
columns_to_clean_2 = ['MBHWISE','MBHWISELO','MBHWISEHI','MBHWISEPUB','MBHWISESOURCE','MBHWISEUPLIM','MSTARWISE','BULGEFRACWISE','MBHWISEQUALITY']

for col in columns_to_clean_1:
    ether[col][to_clean_id] = no_data

for col in columns_to_clean_2:
    if col=='MBHWISEQUALITY':
        ether[col] = ether[col].astype(str)
        ether[col] = str(no_data)[:7]
    else:
        ether[col] = no_data

print('Values in ETHER columns cleaned: {}'.format(len(to_clean_id)))
print('All WISE2MBH columns were cleaned!')
# %%
size_before_depure = len(input_sample)
print('Depuring AllWISE subsample ({})...'.format(directory))
print('Size before depuring {}: {}'.format(directory,size_before_depure))
input_sample = unique(input_sample, keys=['RAJ2000','DEJ2000'], keep='first')

input_sample = input_sample[(input_sample['Z']>0) & (input_sample['Z']!=no_data)]                                                  #Masking of null redshifts
qph_list = input_sample['qph'].tolist()                                                         #Extract quality flags
null_phot = ['X', 'Z']                                                                          #Null and uplim values in W bands
not_null_quality = [not any(letter in qph[:-1] for letter in null_phot) for qph in qph_list]    #Test quality flags

input_sample = input_sample[not_null_quality].copy()                                            #Keep non-null values in quality

size_after_depure = len(input_sample)
rejected = size_before_depure - size_after_depure
print('Size after depure: {}'.format(size_after_depure))

input_sample['INTERNAL_ID'] = np.arange(0,len(input_sample))

for err in ['e_W1mag','e_W2mag','e_W3mag']:
    try:
        fill_err_zero = input_sample[err].data.filled(0)
        input_sample[err] = fill_err_zero
        print('{} band present null entries in {} column. Filled with zeros'.format(err[2:4], err))
    except:
        print('Every {} band entry have 1-sigma error in the {} column.'.format(err[2:4], err))
    

input_sample.add_index('INTERNAL_ID')

rows = np.arange(0,len(input_sample)+chunk,chunk)

''' 
Main algorithm
'''
# %%
print('Running WISE2MBH for catalog in {}:'.format(directory))
for index in range(0,len(rows)-1):
    if index==0:
        start_time = time.time()
    else:
        print('Iteration from row {} to {} out of {} rows ({}%)'.format(rows[index], rows[index+1], rows[-1], np.around(((rows[index]+rows[index+1])/2)/rows[-1]*100, decimals=2)))
    
    allwise = input_sample.iloc['INTERNAL_ID',rows[index]:rows[index+1]]                    #Take slice of indexes
    size_sub_sample = len(allwise)

    allwise['W1-W2_obs'] = allwise['W1mag'] - allwise['W2mag']
    allwise['W2-W3_obs'] = allwise['W2mag'] - allwise['W3mag']

    w1 = wm.array_montecarlo(allwise['W1mag'],allwise['e_W1mag'],n=mc_size)
    w2 = wm.array_montecarlo(allwise['W2mag'],allwise['e_W2mag'],n=mc_size)
    w3 = wm.array_montecarlo(allwise['W3mag'],allwise['e_W3mag'],n=mc_size)
    z = np.array(allwise['Z'])

    w1w2_obs = w1.astype(float)-w2.astype(float)
    w2w3_obs = w2.astype(float)-w3.astype(float)

    first_reject_zone = allwise['W2-W3_obs']<=4.4

    w1 = w1[first_reject_zone]
    w2 = w2[first_reject_zone]
    w3 = w3[first_reject_zone]
    z = z[first_reject_zone]

    w1w2_obs = w1w2_obs[first_reject_zone]
    w2w3_obs = w2w3_obs[first_reject_zone]

    allwise = allwise[first_reject_zone]

    w1w2_kcors = np.zeros(np.shape(w1w2_obs)[0])                                                
    w2w3_kcors = np.zeros(np.shape(w2w3_obs)[0])
    f1complete = np.zeros(np.shape(w1)[0])
    
    '''
    Classifying objects between Galaxies or AGN/QSO
    '''

    object_condition = (allwise['AGNTYPE']==3) | (allwise['AGNTYPE']==4) | (allwise['AGNTYPE']==6) | (allwise['AGNTYPE']==7) | (allwise['AGNTYPE']==8) | (allwise['AGNTYPE']==9) | (allwise['AGNTYPE']==10) | (allwise['AGNTYPE']==22)
    color_condition_1 = (allwise['W1-W2_obs']>0.8) & (allwise['W2-W3_obs']<2.2)
    color_condition_2 = (allwise['W1-W2_obs']>wm.w1w2_treshold_qso(allwise['W2-W3_obs'])) & (allwise['W2-W3_obs']>=2.2) & (allwise['W2-W3_obs']<=4.4)

    optimal_cond = (allwise['Z']<0.5) & ~(object_condition | color_condition_1 | color_condition_2)
    suboptimal_cond = (allwise['Z']>=0.5) & (allwise['Z']<=3) & ~(object_condition | color_condition_1 | color_condition_2)
    nok_cond = (allwise['Z']>3) | object_condition | color_condition_1 | color_condition_2

    optimal_sample = allwise[optimal_cond]
    suboptimal_sample = allwise[suboptimal_cond]                                                                            #Samples for k-correcion                
    nok_sample = allwise[nok_cond]
    
    '''
    Estimating SFR from obs. W3 magnitude and W2-W3 color
    '''

    sfr = wm.w3_to_SFR(w3,w2w3_obs,z[:,None],mc=True,n=mc_size)
    ids_cont = np.where(object_condition & (color_condition_1 | color_condition_2))
    sfr_cont = wm.w3_to_SFR(w3[ids_cont],w2w3_obs[ids_cont],z[:,None][ids_cont],ulirgs=True,mc=True,n=mc_size)

    sfr[ids_cont] = sfr_cont

    allwise['logSFR'] = np.median(sfr, axis=1)          #16, 50 and 84 percentile values saved for final SFR
    allwise['low_logSFR'], allwise['high_logSFR'] = np.percentile(sfr, [16,84], axis=1)
    allwise['SFR_Alert'] = 1*(object_condition & (color_condition_1 | color_condition_2))
    allwise['SFR_Alert'] = np.where((allwise['Z']>0.5),2,allwise['SFR_Alert'])
    
    '''
    Calculating K-corrections for W1, W2 and W3 magnitudes
    '''
                                                                                           
    allwise['K_QUALITY'] = 0
    allwise['K_QUALITY'] = np.where(suboptimal_cond, 1, allwise['K_QUALITY'])
    allwise['K_QUALITY'] = np.where(nok_cond, 2, allwise['K_QUALITY'])

    e_table = wm.kcorr_table('E')
    l_table = wm.kcorr_table('L')
    s_table = wm.kcorr_table('S')

    factor_w2 = 171.787
    factor_w3 = 31.674

    w2w3_limit = -2.5*(np.log10(factor_w3/factor_w2) + 0.1)                             #Last value is the flux limit from Mateos et al. (2012)

    identifier = np.random.choice(np.arange(10), size=3, replace=False)
    for idn, sample, cond in zip(identifier,[optimal_sample, suboptimal_sample, nok_sample], [optimal_cond,suboptimal_cond,nok_cond]):
        if len(sample)!=0:
            if idn!=identifier[2]:                                                                      #Color k-correction for sample with z<0.5
                object_condition = (sample['AGNTYPE']==3) | (sample['AGNTYPE']==4) | (sample['AGNTYPE']==6) | (sample['AGNTYPE']==7) | (sample['AGNTYPE']==8) | (sample['AGNTYPE']==9) | (sample['AGNTYPE']==10) | (sample['AGNTYPE']==22)

                e_index = np.where((sample['T']<=-3) & ~object_condition)[0]
                l_index = np.where((sample['T']>-3) & (sample['T']<=0) & ~object_condition)[0]          #Divide samples in Elliptical, Lenticular and Spiral by T value
                s_index = np.where((sample['T']>0) & (sample['T']!=no_data) & ~object_condition)[0]

                no_type_index = np.where((sample['T']==no_data) & ~object_condition)[0]                 #Sample without T values
                w1w2_no_type = sample[no_type_index]['W1-W2_obs']                                      
                w2w3_no_type = sample[no_type_index]['W2-W3_obs']  
                z_no_type = sample[no_type_index]['Z']                                                    

                e_f1_kcorrected = wm.get_correction_factor(lookup_table=e_table, redshift=sample[e_index]['Z'], correction_factor='f1')
                l_f1_kcorrected = wm.get_correction_factor(lookup_table=l_table, redshift=sample[l_index]['Z'], correction_factor='f1')   #W1 k-corrected 
                s_f1_kcorrected = wm.get_correction_factor(lookup_table=s_table, redshift=sample[s_index]['Z'], correction_factor='f1')

                e_w2w3_kcor = wm.get_correction_factor(lookup_table=e_table, redshift=sample[e_index]['Z'], correction_factor='W2-W3')
                l_w2w3_kcor = wm.get_correction_factor(lookup_table=l_table, redshift=sample[l_index]['Z'], correction_factor='W2-W3')     #W2-W3 k-correction factor
                s_w2w3_kcor = wm.get_correction_factor(lookup_table=s_table, redshift=sample[s_index]['Z'], correction_factor='W2-W3')

                e_w1w2_kcor = wm.get_correction_factor(lookup_table=e_table, redshift=sample[e_index]['Z'], correction_factor='W1-W2')
                l_w1w2_kcor = wm.get_correction_factor(lookup_table=l_table, redshift=sample[l_index]['Z'], correction_factor='W1-W2')     #W1-W2 k-correction factor
                s_w1w2_kcor = wm.get_correction_factor(lookup_table=s_table, redshift=sample[s_index]['Z'], correction_factor='W1-W2')

                no_type_w1w2_kcor = [
                    wm.get_correction_factor(lookup_table=e_table, redshift=z, correction_factor='W1-W2') if (color_x<=w2w3_limit)
                    else wm.get_correction_factor(lookup_table=s_table, redshift=z, correction_factor='W1-W2') if (color_x>w2w3_limit) else 0
                    for color_x,color_y,z in zip(w2w3_no_type,w1w2_no_type, z_no_type)
                    ]
                no_type_w2w3_kcor = [
                    wm.get_correction_factor(lookup_table=e_table, redshift=z, correction_factor='W2-W3') if (color_x<=w2w3_limit)
                    else wm.get_correction_factor(lookup_table=s_table, redshift=z, correction_factor='W2-W3') if (color_x>w2w3_limit) else 0
                    for color_x,color_y,z in zip(w2w3_no_type,w1w2_no_type, z_no_type)
                    ]

                no_type_f1_kcor = [
                    wm.get_correction_factor(lookup_table=e_table, redshift=z, correction_factor='f1') if (color_x<=w2w3_limit)
                    else wm.get_correction_factor(lookup_table=s_table, redshift=z, correction_factor='f1') if (color_x>w2w3_limit) else 1
                    for color_x,color_y,z in zip(w2w3_no_type,w1w2_no_type, z_no_type)
                    ]

                w1w2_kcors[np.where(cond)[0][e_index]] = e_w1w2_kcor
                w1w2_kcors[np.where(cond)[0][l_index]] = l_w1w2_kcor                   #W1-W2 k-correctrion factors filled
                w1w2_kcors[np.where(cond)[0][s_index]] = s_w1w2_kcor
                w1w2_kcors[np.where(cond)[0][no_type_index]] = no_type_w1w2_kcor

                w2w3_kcors[np.where(cond)[0][e_index]] = e_w2w3_kcor
                w2w3_kcors[np.where(cond)[0][l_index]] = l_w2w3_kcor                   #W2-W3 k-correctrion factors filled
                w2w3_kcors[np.where(cond)[0][s_index]] = s_w2w3_kcor
                w2w3_kcors[np.where(cond)[0][no_type_index]] = no_type_w2w3_kcor

                f1complete[np.where(cond)[0][e_index]] = e_f1_kcorrected
                f1complete[np.where(cond)[0][l_index]] = l_f1_kcorrected                 #W1abs k-corrected filled
                f1complete[np.where(cond)[0][s_index]] = s_f1_kcorrected
                f1complete[np.where(cond)[0][no_type_index]] = no_type_f1_kcor
                
            else:
                w1w2_kcors[np.where(cond)[0]] = 0
                w2w3_kcors[np.where(cond)[0]] = 0                                   #For sources with z>3, no k-correction is applied
                f1complete[np.where(cond)[0]] = 1
        else:
            continue
 
    w1w2_kcorrected = w1w2_obs - w1w2_kcors[:,None]                                 #Colors for complete slice are applied
    w2w3_kcorrected = w2w3_obs - w2w3_kcors[:,None]
    
    w1_flux = 309.540 * (10**(-w1/2.5))

    w1_flux_k_corrected = w1_flux*f1complete[:,None]

    w1_ab_mag = -2.5*np.log10(w1_flux_k_corrected) + 8.926
    w1_vega_mag_k_corrected = w1_ab_mag-2.699

    w1abs = w1_vega_mag_k_corrected - wm.distance_modulus_z(z)[:,None]

    allwise['W1-W2_kcor_f'] = w1w2_kcors
    allwise['W2-W3_kcor_f'] = w2w3_kcors

    allwise['W1-W2_kcor'] = np.median(w1w2_kcorrected, axis=1)                      #Median value is saved in data frame
    allwise['W2-W3_kcor'] = np.median(w2w3_kcorrected, axis=1)
    
    '''
    Estimating total stellar mass from corrected W1 magnitude and W1-W2 color
    '''

    w1w2_sat_top = wm.clipping_dist(w1w2_kcorrected, 0.6)                           #W1-W2 is saturated between -0.2 and 0.6 for M/L ratios
    w1w2_sat_complete = wm.clipping_dist(w1w2_sat_top, -0.2, greater_than=False)

    allwise['W1_abs'] = np.median(w1abs, axis=1)
    allwise['W1-W2_clipped'] = np.median(w1w2_sat_complete, axis=1)

    log_sm = wm.wise_to_logsm(w1abs, w1w2_sat_complete)

    allwise['prior_logSM'] = np.median(log_sm, axis=1)

    log_sm_res = wm.wise_to_logsm(w1abs, w1w2_sat_complete, resolved=True)          #log of Stellar mass is calculated and its median is saved

    allwise['prior_logSMres'] = np.median(log_sm_res, axis=1)

    log_sm[np.where(((allwise['ex']==b'5') | (allwise['ex']==b'4')) & (allwise['W2-W3_kcor']<w2w3_limit))[0]] = log_sm_res[np.where(((allwise['ex']==b'5') | (allwise['ex']==b'4')) & (allwise['W2-W3_kcor']<w2w3_limit))[0]]

    allwise['logSM'] = np.median(log_sm, axis=1)
    
    low_sm = 6.5            
    high_sm = 13                                                                    #Limit values in log stellar and masking
                                                                                    
    log_sm_mask = (allwise['logSM']>=low_sm) & (allwise['logSM']<=high_sm)

    log_sm = np.delete(log_sm, np.where(~log_sm_mask)[0], axis=0)                   #logsm and w2w3 gaussian arrays are masked by the limits
    w2w3_kcorrected = np.delete(w2w3_kcorrected, np.where(~log_sm_mask)[0], axis=0)
    allwise = allwise[log_sm_mask]                                                  #data frame is masked

    object_condition = (allwise['AGNTYPE']==3) | (allwise['AGNTYPE']==4) | (allwise['AGNTYPE']==6) | (allwise['AGNTYPE']==7) | (allwise['AGNTYPE']==8) | (allwise['AGNTYPE']==9) | (allwise['AGNTYPE']==10) | (allwise['AGNTYPE']==22)
    color_condition_1 = (allwise['W1-W2_kcor']>0.8) & (allwise['W2-W3_kcor']<2.2)
    color_condition_2 = (allwise['W1-W2_kcor']>wm.w1w2_treshold_qso(allwise['W2-W3_kcor'])) & (allwise['W2-W3_kcor']>=2.2) & (allwise['W2-W3_kcor']<=4.4)
    
    allwise_estim_cond = ~(object_condition | color_condition_1 | color_condition_2)
    allwise_uplim_cond = object_condition | color_condition_1 | color_condition_2
    
    '''
    AGN compensation
    '''

    allwise['BT'] = no_data
    allwise['BT'] = np.where(allwise_uplim_cond, 1, allwise['BT'])

    w1w2_for_agn = allwise['W1-W2_obs'].value.copy()
    w1w2_for_agn[w1w2_for_agn>1.2] = 1.2
    w1w2_for_agn[w1w2_for_agn<0.5] = 0.5

    agn_frac_dirty = wm.agn_fraction(w1w2_for_agn)
    agn_frac_ready = np.where(allwise_estim_cond, 1, agn_frac_dirty)

    allwise['AGN_FRACTION'] = agn_frac_ready
    
    '''
    Estimating T-type from corrected W2-W3 color
    '''

    cond_change_t = (allwise_estim_cond) & ((allwise['T']==no_data) | (allwise['T']>8) | (allwise['T']<-5))     #Condition to change morphological value
    w2w3_to_use = w2w3_kcorrected[np.where(cond_change_t)[0]]                                                                   #w2w3 is masked once again, now for noQSO that require a new T
    new_t_value = wm.w2w3_to_morph(w2w3_to_use)                                                                                 #new T values are calculated using S-shape curve

    t_value_dist = wm.array_montecarlo(allwise['T'], np.zeros(len(allwise['T'])),n=mc_size)  #Gaussians of every T value
    t_value_dist[cond_change_t] = new_t_value

    allwise['T_USED'] = -99                                                                                                  #T values selected are changed
    allwise['T_USED'] = np.where(allwise_estim_cond, np.median(t_value_dist, axis=1), allwise['T_USED'])

    allwise['T_QUALITY'] = 0
    allwise['T_QUALITY'] = np.where(cond_change_t, 1, allwise['T_QUALITY'])
    allwise['T_QUALITY'] = np.where(allwise_uplim_cond, 2, allwise['T_QUALITY'])

    '''
    Obtaining B/T from T-type, obtaining bulge mass from B/T and total stellar mass and then MBH
    '''

    bulge_frac = wm.morph_to_bulge_ratio(t_value_dist)                                               #Bulge fractions are calculated for noQSO using the T value
    allwise['BT'] = np.where(allwise_estim_cond, np.median(bulge_frac, axis=1), allwise['BT'])

    bf_all = wm.array_montecarlo(np.ones(len(allwise)), np.zeros(len(allwise)),n=mc_size)    #Bulge ratios are overwritten and its median is saved
    bf_all[allwise_estim_cond] = bulge_frac[allwise_estim_cond]
    
    log_bm = np.log10(bf_all) + log_sm                                  #Bulge mass is calculated with Stellar mass
    log_sm_agn_cleaned = np.log10(allwise['AGN_FRACTION'])[:,None] + log_bm
  
    log_mbh = wm.bulge_to_mbh(log_bm,mc=True,n=mc_size)
    log_mbh_cleaned = wm.bulge_to_mbh(log_sm_agn_cleaned,mc=True,n=mc_size)

    comp_mbh = wm.comp_mbh(log_mbh)                                      #Empiriclly ompensated MBH
    comp_mbh_cleaned = wm.comp_mbh(log_mbh_cleaned)

    allwise['logMBH'] = np.median(comp_mbh, axis=1)          #16, 50 and 84 percentile values saved for final MBH
    allwise['low_logMBH'], allwise['high_logMBH'] = np.percentile(comp_mbh, [16,84], axis=1)

    allwise['logMBH_AGN_Cleaned'] = np.median(comp_mbh_cleaned, axis=1)          #16, 50 and 84 percentile values saved for final MBH
    allwise['low_logMBH_AGN_Cleaned'], allwise['high_logMBH_AGN_Cleaned'] = np.percentile(comp_mbh_cleaned, [16,84], axis=1)

    allwise['MBHWISEUPLIM'] = 0
    allwise['MBHWISEUPLIM'] = np.where(allwise_uplim_cond, 1, allwise['MBHWISEUPLIM'])

    allwise['MBHWISEQUALITY'] = 4
    allwise['MBHWISEQUALITY'] = np.where(allwise_uplim_cond, 9, allwise['MBHWISEQUALITY'])

    phot_info = [x[:-1] for x in allwise['qph']]
    ex_info = allwise['ex']                             #Photometric quality, extention flag and uplim are combined in one string
    mbh_info = allwise['MBHWISEUPLIM'].astype(str)
    k_quality = allwise['K_QUALITY'].astype(str)
    t_quality = allwise['T_QUALITY'].astype(str)

    quality_flag = []
    for i in np.arange(len(phot_info)):
        str_info = phot_info[i] + ex_info[i] + mbh_info[i] + k_quality[i] + t_quality[i]
        quality_flag.append(str_info)

    allwise['QF'] = quality_flag 

    mbh_non_outliers = allwise['logMBH']>=5
    non_reject_sources = allwise['W2-W3_kcor']<=4.4

    final_allwise = allwise[mbh_non_outliers & non_reject_sources]

    size_after_algorithm = len(final_allwise)
    to_add = size_sub_sample - size_after_algorithm

    rejected += to_add

    cols = ['W1-W2_obs','W2-W3_obs','logSFR','low_logSFR','high_logSFR','SFR_Alert','K_QUALITY','W1-W2_kcor_f','W2-W3_kcor_f','W1-W2_kcor','W2-W3_kcor','W1_abs','W1-W2_clipped','prior_logSM','prior_logSMres','logSM','BT','AGN_FRACTION','T_USED','T_QUALITY','logMBH','low_logMBH','high_logMBH', 'logMBH_AGN_Cleaned', 'low_logMBH_AGN_Cleaned', 'high_logMBH_AGN_Cleaned' ,'MBHWISEUPLIM','MBHWISEQUALITY','QF']

    if not np.isin('K_QUALITY',input_sample.colnames):
        for col in cols:
            input_sample[col] = (np.ones(len(input_sample))*no_data).astype(allwise[col].dtype)
 
    input_sample.loc[final_allwise['INTERNAL_ID']] = final_allwise
        
    if index==0:
        end_time = time.time()
        elapsed_time = end_time - start_time
        print('Estimated total time: {} minutes'.format(np.around((elapsed_time*len(rows)/60), decimals=1)))
        print('Iteration from row {} to {} out of {} rows ({}%)'.format(rows[index], rows[index+1], rows[-1], np.around(((rows[index]+rows[index+1])/2)/rows[-1]*100, decimals=2)))
        
'''
Incorporating results into ETHER
'''

input_sample  = input_sample[(input_sample['logMBH']!=no_data) & (input_sample['logMBH']>=5)]

print('Succesfully finished!')
print('Total rejected sources: {}'.format(rejected))
# %%
print('Overwriting ETHER subsample ({})...'.format(ether_directory))
new_mbh_id = (input_sample['MBH']==no_data) | (input_sample['MBHPUB']==50)

master_id = input_sample['SUBSAMPLE_ID'][new_mbh_id]  
cat_id = new_mbh_id 

ether['MBH'][master_id] = input_sample['logMBH'][cat_id]
ether['MBHLO'][master_id] = input_sample['low_logMBH'][cat_id]
ether['MBHHI'][master_id] = input_sample['high_logMBH'][cat_id]
ether['MBHSOURCE'][master_id] = 10
ether['MBHPUB'][master_id] = 50
ether['MBHQUALITY'][master_id] = input_sample['MBHWISEQUALITY'][cat_id]
ether['MBHLIM'][master_id] = input_sample['MBHWISEUPLIM'][cat_id]

master_id = input_sample['SUBSAMPLE_ID']

ether['MBHWISE'][master_id] = input_sample['logMBH']
ether['MBHWISELO'][master_id] = input_sample['low_logMBH']
ether['MBHWISEHI'][master_id] = input_sample['high_logMBH']
ether['MBHWISE_AGNCLEANED'][master_id] = input_sample['logMBH_AGN_Cleaned']
ether['MBHWISELO_AGNCLEANED'][master_id] = input_sample['low_logMBH_AGN_Cleaned']
ether['MBHWISEHI_AGNCLEANED'][master_id] = input_sample['high_logMBH_AGN_Cleaned']
ether['SFRWISE'][master_id] = input_sample['logSFR']
ether['SFRWISELO'][master_id] = input_sample['low_logSFR']
ether['SFRWISEHI'][master_id] = input_sample['high_logSFR']
ether['SFRWISEALERT'][master_id] = input_sample['SFR_Alert']
ether['MBHWISEPUB'][master_id] = 50
ether['MBHWISESOURCE'][master_id] = 10
ether['MBHWISEUPLIM'][master_id] = input_sample['MBHWISEUPLIM']
ether['MSTARWISE'][master_id] = input_sample['logSM']
ether['BULGEFRACWISE'][master_id] = input_sample['BT']
ether['MBHWISEQUALITY'][master_id] = input_sample['QF']

new_mbh = len(ether[ether['MBHPUB']==50]['MBH'])
total_mbh = len(ether[ether['MBHWISEPUB']==50]['MBHWISE'])

print('New MBH estimations from WISE: {}'.format(new_mbh))
print('Total MBH estimations from WISE: {}'.format(total_mbh))
# %%
ether.write(ether_directory, overwrite=True)
# %%
