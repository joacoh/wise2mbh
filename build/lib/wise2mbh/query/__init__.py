from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack
from astroquery.ipac.ned import Ned
import astropy.units as u
from astroquery.xmatch import XMatch
import numpy as np
import time

sleep_time = 0.2

def query_ned(input_table, ra_column='RA', dec_column='DEC', radius=3, equinox='J2000.0', chunk=400, sleep_time=sleep_time):
    """
    Query NED for additional information based on RA and DEC coordinates or name and merge the results with the input table.

    Inputs:
        table (astropy.table.Table): Input table containing RA and DEC coordinates.
        ra_column (str, optional): Name of the RA column in the table.
        dec_column (str, optional): Name of the DEC column in the table.
        radius (float, optional): Search radius in arcseconds.
        equinox (str, optional): Equinox for coordinates.
        verbose (boolean, optional): If you want to display skiped sources in case of names.
        chunk (int, optional): Slice size to avoid overflow.

    Output:
         table (astropy.table.Table, optinal): Merged table containing NED information. It overwrites the input table.
    """

    table = input_table

    if (ra_column not in table.colnames or dec_column not in table.colnames):
        raise ValueError(f"Columns '{ra_column}' and '{dec_column}' must exist in the input table.")
    
    table['INTERNAL_INDEX'] = list(range(len(table)))
    table['Z'] = np.float64(0)
    table['NED_TYPE'] = 'Unknown'

    if len(table)<=chunk:

        coordinates = SkyCoord(table[ra_column].value*u.deg, table[dec_column].value*u.deg, equinox=equinox)

        for i, coord in enumerate(coordinates):
            ned_results = Ned.query_region(coord, radius=radius*u.arcsec)
            
            if len(ned_results)==0:  
                continue

            ned_results.sort('Separation')
            ned_results['Redshift'] = ned_results['Redshift'].filled(0)
            ned_results = ned_results[0][['Redshift','Type']]
            
            if table['INTERNAL_INDEX'][i] == i:
                table[i][['Z','NED_TYPE']] = ned_results
        
    else:  
        print('Table size bigger than defined chunk ({}), applying slicing.'.format(chunk))
        rows = np.arange(0,len(table)+chunk,chunk)
        ini_query = table[rows[0]:rows[1]]

        coordinates = SkyCoord(ini_query[ra_column].value*u.deg, ini_query[dec_column].value*u.deg, equinox=equinox)

        for i, coord in enumerate(coordinates):
            ned_results = Ned.query_region(coord, radius=radius*u.arcsec)
            
            if len(ned_results)==0:  
                continue

            ned_results.sort('Separation')
            ned_results['Redshift'] = ned_results['Redshift'].filled(0)
            ned_results = ned_results[0][['Redshift','Type']]
            
            if ini_query['INTERNAL_INDEX'][i] == i:
                ini_query[i][['Z','NED_TYPE']] = ned_results       

        for index in range(1, len(rows)-1):
            temp_table = table[rows[index]:rows[index+1]]
            coordinates = SkyCoord(temp_table[ra_column].value*u.deg, temp_table[dec_column].value*u.deg, equinox=equinox)

            for i, coord in enumerate(coordinates):
                ned_results = Ned.query_region(coord, radius=radius*u.arcsec)
                
                if len(ned_results)==0:  
                    continue

                ned_results.sort('Separation')
                ned_results['Redshift'] = ned_results['Redshift'].filled(0)
                ned_results = ned_results[0][['Redshift','Type']]
                
                if temp_table['INTERNAL_INDEX'][i] == i:
                    temp_table[i][['Z','NED_TYPE']] = ned_results      
            
            ini_query = vstack([ini_query,temp_table])
            time.sleep(sleep_time)
        table = ini_query

    table.remove_column('INTERNAL_INDEX')

    object_types_to_keep = ['G', 'QSO', 'RadioS']
    for row in table:
        if row['NED_TYPE'] not in object_types_to_keep:
            row['NED_TYPE'] = 'Unknown'

    return table

def xmatch_by_chunk(input_table, ra_column='RA', dec_column='DEC', radius=3, chunk=400, sleep_time=sleep_time, catalog='vizier:II/328/allwise'):
    """
    XMatch to AllWISE catalog (can be changed) in Vizier to obtain all necesary data for the algorithm based in RA and DEC coordinates.

    Inputs:
        input_table (astropy.table.Table): Input table containing RA and DEC coordinates.
        ra_column (str, optional): Name of the RA column in the table.
        dec_column (str, optional): Name of the DEC column in the table.
        radius (float, optional): Search radius in arcseconds.
        chunk (int, optional): Slice size to avoid overflow.

    Output:
         table (astropy.table.Table, optinal): Merged table containing AllWISE information.
    """

    rows = np.arange(0,len(input_table)+chunk,chunk)
    print('Slicing into {} slices of size {}'.format(len(rows)-1, chunk))
    ini_query = XMatch.query(cat1=input_table[rows[0]:rows[1]],
                    cat2=catalog,
                    max_distance=radius*u.arcsec, colRA1=ra_column,
                    colDec1=dec_column)

    for index in range(1, len(rows)-1):
        print('Iteration from row {} to {} out of {} rows ({}%)'.format(rows[index], rows[index+1], rows[-1], np.around(((rows[index]+rows[index+1])/2)/rows[-1]*100, decimals=2)))
        next_query = XMatch.query(cat1=input_table[rows[index]:rows[index+1]],
                        cat2=catalog,
                        max_distance=radius*u.arcsec, colRA1=ra_column,
                        colDec1=dec_column)
        ini_query = vstack([ini_query,next_query])
        time.sleep(sleep_time)
    
    return ini_query