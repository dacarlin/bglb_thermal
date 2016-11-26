with open('bgl-num.txt') as fn:
    fil = fn.readlines()

import re
counter = 0
for lin in fil:
    spl = lin.split()   
    if re.match(r'[a-z]', spl[0]):
        counter += 1
        spl += [ str(counter) ]
        print(' '.join(spl))
