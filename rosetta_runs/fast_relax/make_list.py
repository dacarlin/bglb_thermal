import pandas 
from Bio.SeqUtils import seq3
df = pandas.read_csv( 'mutant_list', header=None ) 
flag = '-suffix _{} -parser:script_vars target={} new_res={}'
df[0] = df[0].map( lambda x: flag.format( x, x[1:-1], seq3(x[-1]).upper() ) ) 
df.to_csv( 'list', index=False ) 
