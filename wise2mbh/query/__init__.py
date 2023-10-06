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
        coordinates = SkyCoord(table[ra_column]*u.deg, table[dec_column]*u.deg, equinox=equinox)

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
# %%
