'''
WISE2MBH-1.0.1 Sample Creator

This is a general use script to create and format samples to use as input in the pipeline.py script. 

For a tutorial and details, please see the notebook!
'''

import wise2mbh as wm
import numpy as np
from astropy.table import Table

out_dir = 'formated_sample.fits'

#REMOVE THIS 3 LINES OF CODE AND REPLACE THEM WITH YOUR SAMPLE OF NAME, RA AND DEC
test_complete = Table.read('samples/AllWISE_sample.fits')
test_complete = test_complete[np.arange(77512,77712)]
test_to_use = test_complete[['ONAME','RA','DEC']]
#---------------------------------------------------------------------------------

test_ned = wm.query_ned(test_to_use, radius=3)
test_complete = wm.xmatch_by_chunk(test_ned, radius=3)

test_complete.write(out_dir)