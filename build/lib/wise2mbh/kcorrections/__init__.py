import pkg_resources
import pandas as pd

# Load your tables here
ekcorr = pd.read_table(pkg_resources.resource_filename('E.kcorrection_fw1_w1w2_w1w3.tbl'), delim_whitespace=True)
skcorr = pd.read_table(pkg_resources.resource_filename('Sb.kcorrection_fw1_w1w2_w1w3.tbl'), delim_whitespace=True)
lkcorr = pd.read_table(pkg_resources.resource_filename('S0.kcorrection_fw1_w1w2_w1w3.tbl'), delim_whitespace=True)