import fftoolbox as fftb
import numpy as np
from scipy.optimize import curve_fit


import logging, sys


logger = logging.getLogger(__name__)
lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
logging.basicConfig(level=logging.DEBUG, format=lf)

complex_name = 'cf3sh-water'

data = {
        'name': complex_name,
        }

for i in xrange(30):
    parser = fftb.XYZ()
    parser.read_file(filename='data/cf3sh-water/%i.xyz' % i)

    data1, data2 = parser.split_data([1,2,3,4,5,6], [7,8,9])
    data1['name'] = 'cf3sh'
    data2['name'] = 'water'

    data1.update(data)
    data2.update(data)

    molecule1 = fftb.HybridMolecule(data1)

    molecule2 = fftb.HybridMolecule(data2)
    for site in molecule2.sites:
        print site

