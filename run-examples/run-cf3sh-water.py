import numpy as np
from scipy.optimize import curve_fit
import logging, sys
import matplotlib.pyplot as plt
#import seaborn

import fftoolbox as fftb

#logger = logging.getLogger(__name__)
#lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
#logging.basicConfig(level=logging.DEBUG, format=lf)

# get water energy
data = {
        'name': 'water',
        'theory': 'mp2_augccpvtz',
        'representation': ('cartesian', 2),
        'symmetry': False,
        }
log_data = fftb.Gaussian(data=data).data
water_e = log_data['ground energy']

# get CF3SH energy
data['name'] = 'cf3sh'
log_data = fftb.Gaussian(data=data).data
cf3sh_e = log_data['ground energy']

ff = fftb.ForceFieldXML()

DFFEnergies = []
FFEnergies = []
QMEnergies = []

elec = []
lj = []

R = []
d = {
        'C': (0, 0.5),
        'F': (0, 0.5),
        'S': (2, 0.5),
        'H': (0, 0.5),
        }
r = 3.95
for ifile in xrange(30):
    parser = fftb.Gaussian()
    parser.read_file(filename='../data/cf3sh-water/mp2_augccpvtz/%i.log' % ifile, orientation='input')
    total_energy = parser.data['ground energy']
    E = (total_energy - water_e - cf3sh_e)*fftb.au_to_kcal
    QMEnergies.append(E)
    R.append(r)

    cf3sh_data, water_data = parser.split_data([1,2,3,4,5,6], [7,8,9])
    water_data['name'] = 'water'
    water_data['symmetry'] = True

    water = fftb.HybridMolecule(water_data)
    ff.load_forcefields(filename='tip5p.xml', molecule=water)
    
    cf3sh_data['theory'] = 'mp2_augccpvtz'
    cf3sh_data['name'] = 'cf3sh'
    cf3sh_data['symmetry'] = False
    cf3sh_data['sphere params'] = d
    cf3sh_data['representation'] = ('spherical', 2)
    parser = fftb.GDMA(data=cf3sh_data)
    for i, atom in enumerate(cf3sh_data['atoms']):
        atom['multipoles'] = parser.data['atoms'][i]['multipoles']
    cf3sh_data['multipoles'] = parser.data['multipoles']
    cf3sh = fftb.DistributedLebedevMolecule(data=cf3sh_data)
    ff.load_forcefields(molecule=cf3sh, here=True)

    two_molecules = fftb.Complex(cf3sh, water)
    two_molecules.write_mol2(filename='%i.mol2' % ifile, here=True)
    Eff = two_molecules.point_charge_energy() + two_molecules.lennard_jones_energy()
    DFFEnergies.append(Eff)
    elec.append(two_molecules.decomposition['electrostatic'])
    lj.append(two_molecules.decomposition['LJ'])

    r -= 0.05


plt.ylim([-10, 10])
plt.plot(R, DFFEnergies, marker='s', lw=3, label='Distributed FF')
plt.plot(R, QMEnergies, marker='s', lw=3, label='QM')
plt.plot(R, elec, '--', label='electrostatic')
plt.plot(R, lj, '--', label='Lennard Jones')
plt.legend(loc='upper right')
plt.savefig('cf3sh-water.pdf')

