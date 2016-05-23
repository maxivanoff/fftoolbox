import fftoolbox as fftb
import numpy as np

import logging, sys

logger = logging.getLogger(__name__)
lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
logging.basicConfig(filename='./log', filemode='w', level=logging.DEBUG, format=lf)
#logging.basicConfig(level=logging.DEBUG, format=lf)

# MeSNO multipoles distribution
distr = {
        'H': (0, 0.5),
        'C': (0, 0.5),
        'S': (2, 0.5),
        'N': (2, 0.5),
        'O': (2, 0.5),
        }
# CF3SH multipoles distribution
distr = {
        'C': (0, 0.5),
        'S': (1, 0.5),
        'H': (0, 0.5),
        }
# CF3SH multipoles distribution
distr = {
        'C': (0, 0.5),
        'F': (0, 0.5),
        'S': (1, 0.5),
        'H': (0, 0.5),
        }
distr = {
        'O': (0, 0.5),
        'H': (0, 0.5),
        }

data = {
        'name': 'water',
        'theory': 'b3lyp_augccpvdz',
        'representation': ('spherical', 2),
        'sphere params': distr,
        'symmetry': False,
        }

parser = fftb.GDMA(data=data)
data.update(parser.data)
molecule = fftb.DistributedLebedevMolecule(data=data)
results = fftb.Results(data['multipoles'])

ff = fftb.ForceFieldXML()
ff.write_file(molecule=molecule, here=True)
ff.load_forcefields(molecule=molecule, here=True)
print molecule

molecule.write_xyz('%s_%s.xyz' % (data['name'], data['theory']), here=True)
molecule.write_mol2('%s_%s.mol2' % (data['name'], data['theory']), here=True)
molecule.color_charges('%s_%s.pymol' % (data['name'], data['theory']), xyzname='%s_%s.xyz' % (data['name'], data['theory']))

data['density'] = 1.5
parser = fftb.GaussianCube(data=data)
data.update(parser.data.copy())
vdw_grid = fftb.vdwGrid(data)
charges = fftb.LeastSquaresCharges(molecule, vdw_grid)
charges.sites_to_solution()
charges.sites_to_charges()
results.add(charges)
report = fftb.Report('vdw', charges)
print report
print results



