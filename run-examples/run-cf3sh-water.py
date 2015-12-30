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

    d = {
            'C': (0, 0.5),
            'F': (0, 0.5),
            'S': (2, 0.5),
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

    charges = fftb.LeastSquaresCharges(cf3sh, grid)
    charges.sites_to_solution()
    charges.sites_to_charges()
    results = fftb.Results(parser.data['multipoles'])
    results.add(charges)
    print results

    ff.load_forcefields(filename='tip5p.xml', molecule=water)
    two_molecules = fftb.Complex(cf3sh, water)
    ff.load_forcefields(filename='cf3sh_mp2_augccpvtz.xml', molecule=cf3sh, nocharge=True)
    two_molecules.write_mol2(filename='%i.mol2' % ifile, here=True)
    print water
    print cf3sh
    Eff = two_molecules.point_charge_energy() + two_molecules.lennard_jones_energy()
    DFFEnergies.append(Eff)

    r -= 0.05

plt.ylim([-50, 50])
plt.plot(R, DFFEnergies, marker='s', lw=3, label='Distributed FF')
plt.plot(R, FFEnergies, marker='s', lw=3, label='FF')
plt.plot(R, QMEnergies, marker='s', lw=3, label='QM')
#plt.plot(R, elec, label='electrostatic')
#plt.plot(R, lj, label='Lennard Jones')
plt.legend(loc='upper right')
plt.savefig('cf3sh-water.pdf')

