import pandas 
df = pandas.read_csv( 'mutant_list', header=None ) 
flag = '-parser:script_vars target={} new_res={}'
df[0] = df[0].map( lambda x: flag.format( x[1:-1], x[-1] ) ) 
df.to_csv( 'list', index=False ) 
