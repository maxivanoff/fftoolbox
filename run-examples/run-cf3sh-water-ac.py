import fftoolbox as fftb
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import seaborn


import logging, sys


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


# Fit CF3SH charges
data['density'] = 1.5
parser = fftb.GaussianCube(data=data)
data.update(parser.data.copy())
grid = fftb.vdwGrid(data)
molecule = fftb.HM(data)
ls = fftb.LSC(grid=grid, molecule=molecule) 
ls.solve()

results = fftb.Results(log_data['multipoles'])
results.add(ls)
print results


ff = fftb.ForceFieldXML()

DFFEnergies = []
FFEnergies = []
QMEnergies = []

elec = []
lj = []

R = []
r = 3.95
for ifile in xrange(30):
    parser = fftb.Gaussian()
    parser.read_file(filename='../data/cf3sh-water/mp2_augccpvtz/%i.log' % ifile)
    total_energy = parser.data['ground energy']
    E = (total_energy - water_e - cf3sh_e)*fftb.au_to_kcal
    QMEnergies.append(E)
    #r = float(parser.data['comment'].split()[3])
    R.append(r)

    cf3sh_data, water_data = parser.split_data([1,2,3,4,5,6], [7,8,9])
    cf3sh_data['name'] = 'cf3sh'
    cf3sh_data['symmetry'] = False
    water_data['name'] = 'water'
    water_data['symmetry'] = True

    cf3sh = fftb.HybridMolecule(cf3sh_data)
    cf3sh.charges_to_sites(charges=ls.charges)
    ff.load_forcefields(filename='cf3sh_mp2_augccpvtz.xml', molecule=cf3sh)

    water = fftb.HybridMolecule(water_data)
    ff.load_forcefields(filename='tip5p.xml', molecule=water)
    
    two_molecules = fftb.Complex(cf3sh, water)

    print water
    print cf3sh

    Eff = two_molecules.point_charge_energy() + two_molecules.lennard_jones_energy()
    elec.append(two_molecules.decomposition['electrostatic'])
    lj.append(two_molecules.decomposition['LJ'])
    FFEnergies.append(Eff)


    r -= 0.05

plt.ylim([-50, 50])
plt.plot(R, FFEnergies, marker='s', lw=3, label='FF')
plt.plot(R, QMEnergies, marker='s', lw=3, label='QM')
plt.plot(R, elec, '--', label='electrostatic')
plt.plot(R, lj, '--', label='Lennard Jones')
plt.legend(loc='upper right')
plt.savefig('cf3sh-water-ac.pdf')

