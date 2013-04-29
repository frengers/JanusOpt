#!/usr/bin/env python

import numpy as np

taucwepp = np.arange(3, 13, (13-3)/15.0 )
Pmmphr = np.arange(21, 60, (60-21)/15.0 )

for t in taucwepp:
    for p in Pmmphr:
        print 'python objective_function.py --pmmphr %f --taucwepp %f' % (p, t)