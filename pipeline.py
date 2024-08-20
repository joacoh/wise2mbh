'''
WISE2MBH Official Pipeline

This is a general use pipeline for correct estimation of MBH using WISE cataloged data.
This pipeline is intended to work with EVERY SAMPLE GENERATED FROM THE TUTORIALS so please create your input sample from the tutorials in this repository.

Please make sure to un-comment the last line to save your final sample in a real path!
'''
# %%
import wise2mbh as wm 
import numpy as np
from astropy.table import Table
import gc
import time

no_data = 9876543.0
verbose = False
chunk = 400

'''
Import and Cleaning pre-processing
'''

directory = 'samples/AllWISE_sample.fits'
save_csv = False

input_sample = Table.read(directory)

size_before_depure = len(input_sample)
print(f'Size before depure: {size_before_depure}')

allowed_z = input_sample[input_sample['Z']!=0]                                                  #Masking of null redshifts
qph_list = input_sample['qph'].tolist()                                                         #Extract quality flags
null_phot = ['X', 'Z']                                                                          #Null and uplim values in W bands
not_null_quality = [not any(letter in qph[:-1] for letter in null_phot) for qph in qph_list]    #Test quality flags

input_sample = allowed_z[not_null_quality].copy()                                               #Keep non-null values in quality

size_after_depure = len(input_sample)
rejected = size_before_depure - size_after_depure
print(f'Size after depure: {size_after_depure}')

input_sample['INTERNAL_ID'] = np.arange(0,len(input_sample))
input_sample['NED_TYPE'] = 'Unknown'
input_sample['T'] = no_data

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

print('Running WISE2MBH for catalog in {}:'.format(directory))
for index in range(0, len(rows)-1):
    if index==0:
        start_time = time.time()
    else:
        print('Iteration from row {} to {} out of {} rows ({}%)'.format(rows[index], rows[index+1], rows[-1], np.around(((rows[index]+rows[index+1])/2)/rows[-1]*100, decimals=2)))
    
    allwise = input_sample.iloc['INTERNAL_ID',rows[index]:rows[index+1]]                    #Take slice of indexes
    size_sub_sample = len(allwise)

    allwise['W1-W2_obs'] = allwise['W1mag'] - allwise['W2mag']
    allwise['W2-W3_obs'] = allwise['W2mag'] - allwise['W3mag']

    w1 = wm.array_montecarlo(allwise['W1mag'],allwise['e_W1mag'])
    w2 = wm.array_montecarlo(allwise['W2mag'],allwise['e_W2mag'])
    w3 = wm.array_montecarlo(allwise['W3mag'],allwise['e_W3mag'])
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
    w1abs = np.zeros(np.shape(w1))

    object_condition = (allwise['NED_TYPE']=='RadioS') | (allwise['NED_TYPE']=='QSO')
    color_condition_1 = (allwise['W1-W2_obs']>0.8) & (allwise['W2-W3_obs']<2.2)
    color_condition_2 = (allwise['W1-W2_obs']>wm.w1w2_treshold_qso(allwise['W2-W3_obs'])) & (allwise['W2-W3_obs']>=2.2) & (allwise['W2-W3_obs']<=4.4)

    suboptimal_cond = (allwise['Z']>=0.5) & (allwise['Z']<=3) & ~object_condition & ~(color_condition_1 | color_condition_2)
    nok_cond = (allwise['Z']>3) | object_condition | color_condition_1 | color_condition_2

    optimal_sample = allwise[(allwise['Z']<0.5) & ~object_condition & ~(color_condition_1 | color_condition_2)]
    suboptimal_sample = allwise[suboptimal_cond]                                                                            #Samples for k-correcion                
    nok_sample = allwise[nok_cond]
                                                                                           
    allwise['K_QUALITY'] = 0
    allwise['K_QUALITY'] = np.where(suboptimal_cond, 1, allwise['K_QUALITY'])
    allwise['K_QUALITY'] = np.where(nok_cond, 2, allwise['K_QUALITY'])

    e_table = wm.kcorr_table('E')
    l_table = wm.kcorr_table('L')
    s_table = wm.kcorr_table('S')

    identifier = np.random.choice(np.arange(10), size=3, replace=False)
    for idn, sample in zip(identifier,[optimal_sample, suboptimal_sample, nok_sample]):
        if len(sample)!=0:
            if idn!=identifier[2]:                                                                      #Color k-correction for sample with z<0.5
                object_condition = (sample['NED_TYPE']=='RadioS') | (sample['NED_TYPE']=='QSO')

                e_index = np.where((sample['T']<=-3) & ~object_condition)[0]
                l_index = np.where((sample['T']>-3) & (sample['T']<=0) & ~object_condition)[0]          #Divide samples in Elliptical, Lenticular and Spiral by T value
                s_index = np.where((sample['T']>0) & (sample['T']!=no_data) & ~object_condition)[0]

                no_type_index = np.where((sample['T']==no_data) & ~object_condition)[0]                 #Sample without T values
                w1w2_no_type = sample[no_type_index]['W1-W2_obs']                                      
                w2w3_no_type = sample[no_type_index]['W2-W3_obs']  
                z_no_type = sample[no_type_index]['Z']                                                    

                e_w1_kcorrected = wm.w1_k_corrected(lookup_table=e_table,w1=w1[e_index],z=z[e_index])
                l_w1_kcorrected = wm.w1_k_corrected(lookup_table=l_table,w1=w1[l_index],z=z[l_index])   #W1 k-corrected 
                s_w1_kcorrected = wm.w1_k_corrected(lookup_table=s_table,w1=w1[s_index],z=z[s_index])

                e_w2w3_kcor = wm.get_correction_factor(lookup_table=e_table, redshift=sample[e_index]['Z'], correction_factor='W2-W3')
                l_w2w3_kcor = wm.get_correction_factor(lookup_table=l_table, redshift=sample[l_index]['Z'], correction_factor='W2-W3')     #W2-W3 k-correction factor
                s_w2w3_kcor = wm.get_correction_factor(lookup_table=s_table, redshift=sample[s_index]['Z'], correction_factor='W2-W3')

                e_w1w2_kcor = wm.get_correction_factor(lookup_table=e_table, redshift=sample[e_index]['Z'], correction_factor='W1-W2')
                l_w1w2_kcor = wm.get_correction_factor(lookup_table=l_table, redshift=sample[l_index]['Z'], correction_factor='W1-W2')     #W1-W2 k-correction factor
                s_w1w2_kcor = wm.get_correction_factor(lookup_table=s_table, redshift=sample[s_index]['Z'], correction_factor='W1-W2')

                #Both k-correction factors for no type sample, using division by observed W2-W3 color
                factor_w2 = 171.787
                factor_w3 = 31.674

                w2w3_limit = -2.5*(np.log10(factor_w3/factor_w2) + 0.1)                             #Last value is the flux limit from Mateos et al. (2012)

                no_type_w1w2_kcor = [
                    wm.get_correction_factor(lookup_table=e_table, redshift=z, correction_factor='W1-W2') if (color_x<=w2w3_limit)
                    else wm.get_correction_factor(lookup_table=s_table, redshift=z, correction_factor='W1-W2') if (color_x>w2w3_limit) else 0
                    for color_x,color_y,z in zip(w2w3_no_type,w1w2_no_type, z_no_type)
                    ]
                no_type_w2w3_kcor = [
                    wm.get_correction_factor(lookup_table=e_table, redshift=z, correction_factor='W2-W3') if (color_x<=w2w3_limit)
                    else wm.get_correction_factor(lookup_table=s_table, redshift=z, correction_factor='W2-W3') if (color_x>w2w3_limit)  else 0
                    for color_x,color_y,z in zip(w2w3_no_type,w1w2_no_type, z_no_type)
                    ]

                no_type_w1_kcorrected = [
                    wm.w1_k_corrected(lookup_table=e_table,w1=w1[no_type_id],z=np.array([z]))[0] if (color_x<=w2w3_limit)
                    else wm.w1_k_corrected(lookup_table=s_table,w1=w1[no_type_id],z=np.array([z]))[0] if (color_x>w2w3_limit) else w1[no_type_index]
                    for color_x,color_y,z,no_type_id in zip(w2w3_no_type,w1w2_no_type, z_no_type, no_type_index)
                    ]   

                e_w1abs_kcorrected = e_w1_kcorrected - wm.distance_modulus_z(sample[e_index]['Z'])[:,None]          #W1 k-correction for e_sources
                l_w1abs_kcorrected = l_w1_kcorrected - wm.distance_modulus_z(sample[l_index]['Z'])[:,None]          #W1 k-correction for l_sources
                s_w1abs_kcorrected = s_w1_kcorrected - wm.distance_modulus_z(sample[s_index]['Z'])[:,None]          #W1 k-correction for s_sources
                no_type_w1abs_kcorrected = no_type_w1_kcorrected - wm.distance_modulus_z(z_no_type)[:,None]         #W1 k-correction for no_type_sources

                if np.shape(no_type_w1abs_kcorrected)==(0,0):
                    no_type_w1abs_kcorrected = no_type_w1abs_kcorrected.reshape(0,1000)

                w1w2_kcors[e_index] = e_w1w2_kcor
                w1w2_kcors[l_index] = l_w1w2_kcor                   #W1-W2 k-correctrion factors filled
                w1w2_kcors[s_index] = s_w1w2_kcor
                w1w2_kcors[no_type_index] = no_type_w1w2_kcor

                w2w3_kcors[e_index] = e_w2w3_kcor
                w2w3_kcors[l_index] = l_w2w3_kcor                   #W2-W3 k-correctrion factors filled
                w2w3_kcors[s_index] = s_w2w3_kcor
                w2w3_kcors[no_type_index] = no_type_w2w3_kcor

                w1abs[e_index] = e_w1abs_kcorrected
                w1abs[l_index] = l_w1abs_kcorrected                 #W1abs k-corrected filled
                w1abs[s_index] = s_w1abs_kcorrected
                w1abs[no_type_index] = no_type_w1abs_kcorrected
                
            else:
                w1abs_kcorrected = w1[np.where(nok_cond)[0]] - wm.distance_modulus_z(z[np.where(nok_cond)[0]])[:,None]
                w1abs[np.where(nok_cond)[0]] = w1abs_kcorrected

                w1w2_kcors[np.where(nok_cond)[0]] = 0
                w2w3_kcors[np.where(nok_cond)[0]] = 0                               #For sources with z>3, no k-correction is applied
        else:
            continue
 
    w1w2_kcorrected = w1w2_obs - w1w2_kcors[:,None]                                 #Colors for complete slice are applied
    w2w3_kcorrected = w2w3_obs - w2w3_kcors[:,None]

    allwise['W1-W2_kcor'] = np.median(w1w2_kcorrected, axis=1)                      #Median value is saved in data frame
    allwise['W2-W3_kcor'] = np.median(w2w3_kcorrected, axis=1)
    w1w2_sat_top = wm.clipping_dist(w1w2_kcorrected, 0.6)                           #W1-W2 is saturated between -0.2 and 0.6 for M/L ratios
    w1w2_sat_complete = wm.clipping_dist(w1w2_sat_top, -0.2, greater_than=False)

    log_sm = wm.wise_to_logsm(w1abs, w1w2_sat_complete)
    log_sm_res = wm.wise_to_logsm(w1abs, w1w2_sat_complete, resolved=True)                             #log of Stellar mass is calculated and its median is saved

    log_sm[np.where(((allwise['ex']==b'5') | (allwise['ex']==b'4')) & (allwise['W2-W3_kcor']<w2w3_limit))[0]] = log_sm_res[np.where(((allwise['ex']==b'5') | (allwise['ex']==b'4')) & (allwise['W2-W3_kcor']<w2w3_limit))[0]]

    allwise['logSM'] = np.median(log_sm, axis=1)
    
    low_sm = 6.5            
    high_sm = 13                                                                    #Limit values in log stellar and masking
                                                                                    
    log_sm_mask = (allwise['logSM']>=low_sm) & (allwise['logSM']<=high_sm)

    log_sm = np.delete(log_sm, np.where(~log_sm_mask)[0], axis=0)       #logsm and w2w3 gaussian arrays are masked by the limits
    w2w3_kcorrected = np.delete(w2w3_kcorrected, np.where(~log_sm_mask)[0], axis=0)
    allwise = allwise[log_sm_mask]                                                      #data frame is masked

    object_condition = (allwise['NED_TYPE']==b'RadioS') | (allwise['NED_TYPE']==b'QSO')
    color_condition_1 = (allwise['W1-W2_kcor']>0.8) & (allwise['W2-W3_kcor']<2.2)
    color_condition_2 = (allwise['W1-W2_kcor']>wm.w1w2_treshold_qso(allwise['W2-W3_kcor'])) & (allwise['W2-W3_kcor']>=2.2) & (allwise['W2-W3_kcor']<=4.4)
    
    allwise_estim_cond = ~(object_condition | color_condition_1 | color_condition_2)
    allwise_uplim_cond = object_condition | color_condition_1 | color_condition_2

    allwise['BT'] = no_data
    allwise['BT'] = np.where(allwise_uplim_cond,1, allwise['BT'])

    cond_change_t = (allwise_estim_cond) & ((allwise['T']==no_data) | (allwise['T']>8) | (allwise['T']<-5))     #Condition to change morphological value
    w2w3_to_use = w2w3_kcorrected[np.where(cond_change_t)[0]]                                                                   #w2w3 is masked once again, now for noQSO that require a new T
    new_t_value = wm.w2w3_to_morph(w2w3_to_use)                                                                                 #new T values are calculated using S-shape curve

    t_value_dist = wm.array_montecarlo(allwise['T'], np.zeros(len(allwise['T'])))  #Gaussians of every T value
    t_value_dist[cond_change_t] = new_t_value

    allwise['T_USED'] = -99                                                                                                  #T values selected are changed
    allwise['T_USED'] = np.where(allwise_estim_cond, np.median(t_value_dist, axis=1), allwise['T_USED'])

    allwise['T_QUALITY'] = 0
    allwise['T_QUALITY'] = np.where(cond_change_t, 1, allwise['T_QUALITY'])
    allwise['T_QUALITY'] = np.where(allwise_uplim_cond, 2, allwise['T_QUALITY'])

    var_names = [
        'w1', 'w2', 'w3', 'w1abs', 'w1w2', 'w2w3',
        'e_sources', 'l_sources', 's_sources', 'no_type_sources',
        'w1w2_kcorrected', 'w1w2_kcors', 'w1w2_no_type', 'w1w2_sat_complete', 'w1w2_sat_top',
        'w2w3_kcorrected', 'w2w3_kcors', 'w2w3_no_type', 'w2w3_to_use'
    ]
    
    for var_name in var_names:
        try:
            del globals()[var_name]  # Try to delete the variable
            if verbose:
                print(f"Deleted {var_name}")
        except KeyError:
            if verbose:
                print(f"{var_name} doesn't exist, skipping")

    bulge_frac = wm.morph_to_bulge_ratio(t_value_dist)                                               #Bulge fractions are calculated for noQSO using the T value
    allwise['BT'] = np.where(allwise_estim_cond, np.median(bulge_frac, axis=1), allwise['BT'])

    bf_all = wm.array_montecarlo(np.ones(len(allwise)), np.zeros(len(allwise)))    #Bulge ratios are overwritten and its median is saved
    bf_all[allwise_estim_cond] = bulge_frac[allwise_estim_cond]

    log_bm = np.log10(bf_all) + log_sm                                  #Bulge mass is calculated with Stellar mass
  
    log_mbh = wm.bulge_to_mbh(log_bm)

    comp_mbh = wm.comp_mbh(log_mbh)                                      #Empiriclly ompensated MBH 

    allwise['logMBH'] = np.median(comp_mbh, axis=1)          #16, 50 and 84 percentile values saved for final MBH
    allwise['low_logMBH'], allwise['high_logMBH'] = np.percentile(comp_mbh, [16,84], axis=1)[0], np.percentile(comp_mbh, [16,84], axis=1)[-1]

    allwise['MBHWISEUPLIM'] = 0
    allwise['MBHWISEUPLIM'] = np.where(allwise_uplim_cond, 1, allwise['MBHWISEUPLIM'])

    allwise['MBHWISEUPLIM'] = 0
    allwise['MBHWISEUPLIM'] = np.where(allwise_uplim_cond, 1, allwise['MBHWISEUPLIM'])
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

    size_after_algorithm = len(allwise)
    to_add = size_sub_sample - size_after_algorithm

    rejected += to_add

    cols = ['W1-W2_obs','W2-W3_obs','K_QUALITY','W1-W2_kcor','W2-W3_kcor','logSM','BT','T_USED','T_QUALITY','logMBH','low_logMBH','high_logMBH','MBHWISEUPLIM','QF']

    if not np.isin('K_QUALITY',input_sample.colnames):
        for col in cols:
            input_sample[col] = (np.ones(len(input_sample))*no_data).astype(allwise[col].dtype)

    input_sample.loc[final_allwise['INTERNAL_ID']] = final_allwise
    
    del allwise, bf_all, bulge_frac, log_bm, log_mbh, log_sm, comp_mbh, new_t_value, t_value_dist
    gc.collect()
    
    if index==0:
        end_time = time.time()
        elapsed_time = end_time - start_time
        print('Estimated total time: {} minutes'.format(np.around((elapsed_time*len(rows)/60)+3.5, decimals=1)))
        print('Iteration from row {} to {} out of {} rows ({}%)'.format(rows[index], rows[index+1], rows[-1], np.around(((rows[index]+rows[index+1])/2)/rows[-1]*100, decimals=2)))

    if index==range(0, len(rows)-1)[-1]:
        print('Succesfully finished!')
        print('Total rejected sources: {}'.format(rejected))

if save_csv:
    input_sample.write(directory[:-5]+'-w2m.csv')
else:
    input_sample.write(directory[:-5]+'-w2m.fits')
# %%

