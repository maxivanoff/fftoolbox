import numpy as np
from scipy.optimize import curve_fit
import logging, sys
import matplotlib.pyplot as plt
import seaborn

import fftoolbox as fftb

#logger = logging.getLogger(__name__)
#lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
#logging.basicConfig(level=logging.DEBUG, format=lf)
def get_bond(a, b):
    return np.linalg.norm(a - b)

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
data['name'] = 'cis-mesno'
log_data = fftb.Gaussian(data=data).data
mesno_e = log_data['ground energy']

ff = fftb.ForceFieldXML()

DFFEnergies = []
FFEnergies = []
QMEnergies = []

elec = []
lj = []
sn_k = []
no_k = []

R = []
l=2
distr = {
        'H': (0, 0.5),
        'C': (0, 0.5),
        'S': (l, 0.5),
        'N': (l, 0.5),
        'O': (l, 0.5),
        }
r = 3.95
for ifile in xrange(353):
    parser = fftb.XYZ()
    parser.read_file(filename='../../data/mesno-ammonium/ordered/cis-%i.xyz' % ifile)
    total_energy = float(parser.data['comment'].split()[1])
    E = (total_energy - ammonium_e - mesno_e)*fftb.au_to_kcal
    QMEnergies.append(E)
    R.append(r)

    mesno_data, ammonium_data = parser.split_data([1,2,3,4,5,6,7], [8,9,10,11,12])
    ammonium_data['name'] = 'ammonium'
    ammonium_data['theory'] = 'pbe_def2svp'
    ammonium_data['symmetry'] = True

    ammonium = fftb.HybridMolecule(ammonium_data)
    ff.load_forcefields(molecule=ammonium, here=True)
    
    mesno_data['theory'] = 'pbe_def2svp'
    mesno_data['name'] = 'cis-mesno'
    mesno_data['symmetry'] = False
    mesno_data['sphere params'] = distr
    mesno_data['representation'] = ('spherical', 2)
    parser = fftb.GDMA(data=mesno_data)
    for i, atom in enumerate(mesno_data['atoms']):
        atom['multipoles'] = parser.data['atoms'][i]['multipoles']
    mesno_data['multipoles'] = parser.data['multipoles']
    mesno = fftb.DistributedLebedevMolecule(data=mesno_data)
    ff.load_forcefields(molecule=mesno, filename='cis-mesno-ff-l-%i.xml' % l, here=True)

    # internal energy
    coordinates = {}
    for site in mesno.atoms:
        coordinates[site.element] = site.coordinates
    s = coordinates['S']*fftb.units.au_to_angst
    n = coordinates['N']*fftb.units.au_to_angst
    o = coordinates['O']*fftb.units.au_to_angst
    sn = get_bond(s, n)
    no = get_bond(n, o)
    e_no = 733.7*(no - 1.1816)**2
    e_sn = 80.5*(sn - 1.7798)**2
    Eff = e_no + e_sn
    Eff = 0.

    no_k.append(e_no)
    sn_k.append(e_sn)
    
    #print ammonium
    #print mesno

    two_molecules = fftb.Complex(mesno, ammonium)
    #two_molecules.write_mol2(filename='%i.mol2' % ifile, here=True)
    Eff += two_molecules.point_charge_energy() + two_molecules.lennard_jones_energy()
    DFFEnergies.append(Eff)
    elec.append(two_molecules.decomposition['electrostatic'])
    lj.append(two_molecules.decomposition['LJ'])

    r -= 0.05


plt.ylim([-20, 5])
plt.plot(R, DFFEnergies, marker='s', lw=3, label='Distributed FF')
plt.plot(R, QMEnergies, marker='s', lw=3, label='QM')
plt.plot(R, elec, '--', label='electrostatic')
plt.plot(R, lj, '--', label='Lennard Jones')
plt.plot(R, no_k, '--', label='NO bond')
plt.plot(R, sn_k, '--', label='SN bond')
plt.legend(loc='lower right')
plt.savefig('cis-mesno-ammonium-fitted-l-%i.pdf' % l)

