import numpy as np
from scipy.optimize import curve_fit
import logging, sys
import matplotlib.pyplot as plt
import seaborn

import fftoolbox as fftb

#logger = logging.getLogger(__name__)
#lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
#logging.basicConfig(level=logging.DEBUG, format=lf)

# get ammonium energy
data = {
        'name': 'ammonium',
        'theory': 'pbe_def2svp',
        'representation': ('cartesian', 2),
        'symmetry': True,
        }
log_data = fftb.Gaussian(data=data).data
ammonium_e = log_data['ground energy']

# get mesno energy
data['name'] = 'trans-mesno'
log_data = fftb.Gaussian(data=data).data
mesno_e = log_data['ground energy']

ff = fftb.ForceFieldXML()

DFFEnergies = []
FFEnergies = []
QMEnergies = []

elec = []
lj = []

R = []
distr = {
        'H': (0, 0.5),
        'C': (0, 0.5),
        'S': (2, 0.5),
        'N': (2, 0.5),
        'O': (2, 0.5),
        }
r = 3.95
for ifile in xrange(335):
    parser = fftb.XYZ()
    parser.read_file(filename='../data/mesno-ammonium/trans-%i.xyz' % ifile)
    total_energy = float(parser.data['comment'].split()[1])
    E = (total_energy - ammonium_e - mesno_e)*fftb.au_to_kcal
    QMEnergies.append(E)
    R.append(r)

    mesno_data, ammonium_data = parser.split_data([1,2,3,4,9,10,11], [5,6,7,8,12])
    ammonium_data['name'] = 'ammonium'
    ammonium_data['theory'] = 'pbe_def2svp'
    ammonium_data['symmetry'] = True

    ammonium = fftb.HybridMolecule(ammonium_data)
    ff.load_forcefields(molecule=ammonium, here=True)
    
    mesno_data['theory'] = 'pbe_def2svp'
    mesno_data['name'] = 'trans-mesno'
    mesno_data['symmetry'] = False
    mesno_data['sphere params'] = distr
    mesno_data['representation'] = ('spherical', 2)
    parser = fftb.GDMA(data=mesno_data)
    for i, atom in enumerate(mesno_data['atoms']):
        atom['multipoles'] = parser.data['atoms'][i]['multipoles']
    mesno_data['multipoles'] = parser.data['multipoles']
    mesno = fftb.DistributedLebedevMolecule(data=mesno_data)
    ff.load_forcefields(molecule=mesno, here=True)
    
    #print ammonium
    #print mesno

    two_molecules = fftb.Complex(mesno, ammonium)
    #two_molecules.write_mol2(filename='%i.mol2' % ifile, here=True)
    Eff = two_molecules.point_charge_energy() + two_molecules.lennard_jones_energy()
    DFFEnergies.append(Eff)
    elec.append(two_molecules.decomposition['electrostatic'])
    lj.append(two_molecules.decomposition['LJ'])

    r -= 0.05


plt.ylim([-20, 5])
plt.plot(R, DFFEnergies, marker='s', lw=3, label='Distributed FF')
plt.plot(R, QMEnergies, marker='s', lw=3, label='QM')
plt.plot(R, elec, '--', label='electrostatic')
plt.plot(R, lj, '--', label='Lennard Jones')
plt.legend(loc='lower right')
plt.savefig('trans-mesno-ammonium.pdf')

