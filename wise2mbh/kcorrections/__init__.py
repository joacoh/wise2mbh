import pkg_resources
import pandas as pd

# Load your tables here
ekcorr = pd.read_table(pkg_resources.resource_filename('wise2mbh.kcorrections','E.kcorrection_fw1_w1w2_w1w3.tbl'), sep='\s+')
skcorr = pd.read_table(pkg_resources.resource_filename('wise2mbh.kcorrections','Sb.kcorrection_fw1_w1w2_w1w3.tbl'), sep='\s+')
lkcorr = pd.read_table(pkg_resources.resource_filename('wise2mbh.kcorrections','S0.kcorrection_fw1_w1w2_w1w3.tbl'), sep='\s+')