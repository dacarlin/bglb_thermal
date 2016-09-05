import argparse
import re 

parser = argparse.ArgumentParser()
parser.add_argument('msa', help='A global multiple sequence alignment. The first sequence will be used as the reference sequence and all other sequences in the file will be trimmed to match')
args = parser.parse_args()

with open( args.msa ) as msa:
  msa = msa.readlines()

numbers = [ i for i, j in enumerate( list ( msa[1] ) ) if re.match( r'[a-z]', j ) ]

for line in msa:
  if line.startswith(r'>'):
    print( line.strip() )
  else:
    alignment = [ j for i, j in enumerate( list( line ) ) if i in numbers ] 
    print( ''.join( alignment ) ) 


