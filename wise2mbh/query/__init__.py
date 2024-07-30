# %%
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astroquery.ipac.ned import Ned
import astropy.units as u
from astroquery.xmatch import XMatch
import numpy as np

def query_ned(table, ra_column='RA', dec_column='DEC', radius=3, equinox='J2000.0', coords=True, name_column='ONAME', verbose=False):
    """
    Query NED for additional information based on RA and DEC coordinates or name and merge the results with the input table.

    Inputs:
        table (astropy.table.Table): Input table containing RA and DEC coordinates.
        ra_column (str, optional): Name of the RA column in the table.
        dec_column (str, optional): Name of the DEC column in the table.
        radius (float, optional): Search radius in arcseconds.
        equinox (str, optional): Equinox for coordinates.
        coords (boolean, optinal): If you want to use coordinates (True) or names (False).
        name_column (str, optional): Name of the NAMES column in the table.
        verbose (boolean, optional): If you want to display skiped sources in case of names.

    Output:
         table (astropy.table.Table, optinal): Merged table containing NED information. It overwrites the input table.
    """

    if (coords) and (ra_column not in table.colnames or dec_column not in table.colnames):
        raise ValueError(f"Columns '{ra_column}' and '{dec_column}' must exist in the input table.")
    
    if (not coords) and (ra_column not in table.colnames or dec_column not in table.colnames):
        raise ValueError(f"Column '{name_column}' must exist in the input table.")

    table['INTERNAL_INDEX'] = list(range(len(table)))
    table['Z'] = np.float64(0)
    table['NED_TYPE'] = 'Unknown'

    if coords:
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
        for i, name in enumerate(table[name_column]):
            try:
                ned_results = Ned.query_region(name, radius=radius*u.arcsec)
            except:
                if verbose:
                    print(f'No source with name: {name}')
                continue
            
            if len(ned_results)==0:  
                continue

            ned_results.sort('Separation')
            ned_results['Redshift'] = ned_results['Redshift'].filled(0)
            ned_results = ned_results[0][['Redshift','Type']]
            
            if table['INTERNAL_INDEX'][i] == i:
                table[i][['Z','NED_TYPE']] = ned_results
            
    table.remove_column('INTERNAL_INDEX')

    object_types_to_keep = ['G', 'QSO', 'RadioS']
    for row in table:
        if row['NED_TYPE'] not in object_types_to_keep:
            row['NED_TYPE'] = 'Unknown'

    return table

def xmatch_allwise(input_table, ra_column='RA', dec_column='DEC', radius=3):
    """
    XMatch to AllWISE catalog in Vizier to obtain all necesary data for the algorithm based in RA and DEC coordinates.

    Inputs:
        input_table (astropy.table.Table): Input table containing RA and DEC coordinates.
        ra_column (str, optional): Name of the RA column in the table.
        dec_column (str, optional): Name of the DEC column in the table.
        radius (float, optional): Search radius in arcseconds.

    Output:
         table (astropy.table.Table, optinal): Merged table containing AllWISE information.
    """
    table = XMatch.query(cat1=input_table,
                     cat2='vizier:II/328/allwise',
                     max_distance=radius*u.arcsec, colRA1=ra_column,
                     colDec1=dec_column)
    return table

"""
Simple XMatch function between Pandas DataFrames
Inputs:
cat: Catalogue to match
master: Master catalogue to match (put here the biggest of the two cat to match)
names: Name column of the sources in cat
cat_ra: RA column in cat
cat_dec: DEC column in cat
master_ra: RA column in master
master_dec: DEC column in master
tol_coord: Radius limit to match a source
sources: If True, print every matched sources. False by default

Outputs:
cat_ids: List of the rows IDs that match in cat
master_ids: List of the rows IDs that match in master
"""
def simple_xmatch(cat,master,names='NAME',cat_ra='RA',cat_dec='DEC',master_ra='RA',master_dec='DEC', tol_coord=5, sources=False):

    ra = cat[cat_ra]
    dec = cat[cat_dec]
    names = cat[names]

    tol_coord = tol_coord/3600
    cat_ids = []
    master_ids = []
    if len(names)!=len(ra)!=len(dec):
        raise ValueError('Names and coordinates must have the same length')
    else:
        print('Crossmatching with a tolerance of {} arc-second(s) difference'.format(tol_coord*3600))
        count = 0  
        for i in np.arange(len(names)):
            matches = 0
            name = names[i]
            ra2match = ra[i]
            dec2match = dec[i]

            master_masked = master[master[master_dec].between(dec2match-tol_coord,dec2match+tol_coord)]
            ra4cat = master_masked[master_ra]
            dec4cat = master_masked[master_dec]

            dist = np.sqrt((ra2match-ra4cat)**2 + (dec2match-dec4cat)**2)   #Substraction to analyze

            cond1 = (dist<=tol_coord)*(min(dist, default=False)==dist)       #Boolean array of condition after match
            cond2 = [j for j, x in enumerate(cond1) if x]                    #Takes the index of True values only
            if len(cond2)!=0:
                cond2 = [cond2[0]] 
                
            for k in cond2:
                master_ids.append(master_masked.index[k])                  #Save index of match in ETHER
                matches += 1
                count += 1

            if matches>0:
                cat_ids.append(cat.index[i])
                if sources:
                    print('Source {} has a match!'.format(name))
            
            if len(cat_ids)!=len(master_ids):
                raise ValueError('Corrupted Series in DataFrame')

        print('Total matches: {}'.format(count))
        return cat_ids, master_ids
# %%
