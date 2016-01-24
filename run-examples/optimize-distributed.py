import fftoolbox as fftb
import numpy as np

import logging, sys

logger = logging.getLogger(__name__)
lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
logging.basicConfig(filename='./log', filemode='w', level=logging.DEBUG, format=lf)
#logging.basicConfig(level=logging.DEBUG, format=lf)

# CF3SH multipoles distribution
distr = {
        'C': (0, 0.5),
        'F': (0, 0.5),
        'S': (2, 0.5),
        'H': (0, 0.5),
        }
# MeSNO multipoles distribution
distr = {
        'H': (0, 0.5),
        'C': (0, 0.5),
        'S': (2, 0.5),
        'N': (2, 0.5),
        'O': (2, 0.5),
        }

data = {
        'name': 'trans-mesno',
        'theory': 'pbe_def2svp',
        'representation': ('spherical', 2),
        'sphere params': distr,
        'symmetry': False,
        }

parser = fftb.GDMA(data=data)
data.update(parser.data)
molecule = fftb.DistributedLebedevMolecule(data=data)
results = fftb.Results(data['multipoles'])

data['density'] = 1.5
parser = fftb.GaussianCube(data=data)
data.update(parser.data.copy())
vdw_grid = fftb.vdwGrid(data)
charges = fftb.LeastSquaresCharges(molecule, vdw_grid)
charges.sites_to_solution()
charges.sites_to_charges()

def set_x(x):
    i = 0
    for atom in molecule.atoms:
        for l in xrange(atom.rank+1):
            for m in xrange(-l,l+1):
                if m < 0:
                    mname = '%i%is' % (l, abs(m))
                if m == 0:
                    mname = '%i0' % l
                if m > 0:
                    mname = '%i%ic' % (l, m)
                if mname in atom.reference_multipoles.keys():
                    atom.reference_multipoles[mname] = x[i]
                    i += 1
        atom.set_charges()

def rmsd_function(x):
    set_x(x)
    charges.sites_to_solution()
    return charges.rmsd

from scipy.optimize import minimize

x0 = []
for atom in molecule.atoms:
    for l in xrange(atom.rank+1):
        for m in xrange(-l,l+1):
            if m < 0:
                mname = '%i%is' % (l, abs(m))
            if m == 0:
                mname = '%i0' % l
            if m > 0:
                mname = '%i%ic' % (l, m)
            if mname in atom.reference_multipoles.keys():
                x0.append(atom.reference_multipoles[mname])

res = minimize(rmsd_function, x0, method='BFGS')

print res

set_x(res.x)
charges.sites_to_solution()
results.add(charges)
print results

