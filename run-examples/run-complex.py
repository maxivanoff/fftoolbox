import fftoolbox as fftb
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import seaborn


import logging, sys


#logger = logging.getLogger(__name__)
#lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
#logging.basicConfig(level=logging.DEBUG, format=lf)

ff = fftb.ForceFieldXML()
parser = fftb.Gaussian()
parser.read_file(filename='../data/cf3sh-water/mp2_augccpvtz/0.log', orientation='input')

cf3sh_data, water_data = parser.split_data([1,2,3,4,5,6], [7,8,9])
cf3sh_data['name'] = 'cf3sh'
cf3sh_data['symmetry'] = False
water_data['name'] = 'water'
water_data['symmetry'] = True

water = fftb.HybridMolecule(water_data)
ff.load_forcefields(filename='tip5p.xml', molecule=water)

d = {
        'C': (0, 0.5),
        'F': (0, 0.5),
        'S': (1, 0.5),
        'H': (0, 0.5),
        }
cf3sh_data['theory'] = 'mp2_augccpvtz'
cf3sh_data['symmetry'] = False
cf3sh_data['sphere params'] = d
cf3sh_data['representation'] = ('spherical', 2)
parser = fftb.GDMA(data=cf3sh_data)
for i, atom in enumerate(cf3sh_data['atoms']):
    atom['multipoles'] = parser.data['atoms'][i]['multipoles']
cf3sh_data['multipoles'] = parser.data['multipoles']

cf3sh = fftb.DistributedLebedevMolecule(data=cf3sh_data)

q1 = cf3sh.multipoles['11c']
q2 = cf3sh.multipoles['11s']
q3 = cf3sh.multipoles['10']
#print q1,cf3sh_data['multipoles']['11c']
print 'my molecular |Q1|=', np.sqrt(q1*q1 + q2*q2 + q3*q3)
#print cf3sh.frames[0].center.multipoles
q1 = cf3sh_data['multipoles']['11c']
q2 = cf3sh_data['multipoles']['11s']
q3 = 0.0
print 'reference molecular |Q1|=', np.sqrt(q1*q1 + q2*q2 + q3*q3)
q1 = cf3sh.frames[0].center.multipoles['11c']
q2 = cf3sh.frames[0].center.multipoles['11s']
q3 = cf3sh.frames[0].center.multipoles['10']
print 'my S |Q1|=', np.sqrt(q1*q1 + q2*q2 + q3*q3)
q1 = cf3sh_data['atoms'][1]['multipoles']['11c']
q2 = cf3sh_data['atoms'][1]['multipoles']['11s']
q3 = 0.0
print 'reference S |Q1|=', np.sqrt(q1*q1 + q2*q2 + q3*q3)


two_molecules = fftb.Complex(cf3sh, water)
ff.load_forcefields(filename='cf3sh_mp2_augccpvtz.xml', molecule=cf3sh, nocharge=True)
two_molecules.write_mol2(filename='0.mol2', here=True)
print water
print cf3sh
Eff = two_molecules.point_charge_energy() + two_molecules.lennard_jones_energy()
print 'energy', Eff

