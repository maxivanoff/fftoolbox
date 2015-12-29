import fftoolbox as fftb
import numpy as np
from scipy.optimize import curve_fit


import logging, sys


logger = logging.getLogger(__name__)
lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
logging.basicConfig(level=logging.DEBUG, format=lf)

ff = fftb.ForceFieldXML()
data = {
        'name': 'cf3sh',
        'theory': 'mp2_augccpvtz',
        'density': 1.5,
        'representation': ('cartesian', 2)
        }
results = fftb.Results(fftb.Gaussian(data=data).data['multipoles'])
parser = fftb.GaussianCube(data=data)
data.update(parser.data.copy())
grid = fftb.vdwGrid(data)
molecule = fftb.HM(data)
ls = fftb.LSC(grid=grid, molecule=molecule) 
ls.solve()
results.add(ls)
print results

for i in xrange(30):
    parser = fftb.XYZ()
    parser.read_file(filename='../data/cf3sh-water/%i.xyz' % i)

    cf3sh_data, water_data = parser.split_data([1,2,3,4,5,6], [7,8,9])
    cf3sh_data['name'] = 'cf3sh'
    water_data['name'] = 'water'

    cf3sh = fftb.HybridMolecule(cf3sh_data)
    water = fftb.HybridMolecule(water_data)
    
    two_molecules = fftb.Complex(cf3sh, water)

    ff.load_forcefields(filename='tip5p.xml', molecule=water)
    cf3sh.charges_to_sites(charges=ls.charges)

    two_molecules.write_mol2(filename='%i.mol2' % i, here=True)

    print water
    break

