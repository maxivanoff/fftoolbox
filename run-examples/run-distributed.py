import fftoolbox as fftb
import numpy as np

import logging, sys

logger = logging.getLogger(__name__)
lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
logging.basicConfig(filename='./log', filemode='w', level=logging.DEBUG, format=lf)
#logging.basicConfig(level=logging.DEBUG, format=lf)

d = {
        'H': (0, 0.5),
        'C': (0, 0.5),
        'S': (2, 0.5),
        'N': (2, 0.5),
        'O': (2, 0.5),
        }

data = {
        'name': 'cis-mesno',
        'theory': 'mp2_augccpvtz',
        'representation': ('spherical', 2),
        'sphere params': d,
        }

parser = fftb.GDMA(data=data)
data.update(parser.data)
molecule = fftb.DistributedLebedevMolecule(data=data)
results = fftb.Results(data['multipoles'])

#ff = fftb.ForceFieldXML()
#ff.write_file(molecule=molecule)

molecule.write_xyz('%s_%s.xyz' % (data['name'], data['theory']), here=True)
molecule.write_mol2('%s_%s.mol2' % (data['name'], data['theory']), here=True)
molecule.color_charges('%s_%s.pymol' % (data['name'], data['theory']), xyzname='%s_%s.xyz' % (data['name'], data['theory']))


q1 = molecule.multipoles['11c']
q2 = molecule.multipoles['11s']
q3 = molecule.multipoles['10']
#print q1,cf3sh_data['multipoles']['11c']
print 'my molecular |Q1|=', np.sqrt(q1*q1 + q2*q2 + q3*q3)
#print cf3sh.frames[0].center.multipoles
q1 = data['multipoles']['11c']
q2 = data['multipoles']['11s']
q3 = 0.0
print 'reference molecular |Q1|=', np.sqrt(q1*q1 + q2*q2 + q3*q3)
q1 = molecule.frames[0].center.multipoles['11c']
q2 = molecule.frames[0].center.multipoles['11s']
q3 = molecule.frames[0].center.multipoles['10']
print 'my S |Q1|=', np.sqrt(q1*q1 + q2*q2 + q3*q3)
q1 = data['atoms'][1]['multipoles']['11c']
q2 = data['atoms'][1]['multipoles']['11s']
q3 = 0.0
print 'reference S |Q1|=', np.sqrt(q1*q1 + q2*q2 + q3*q3)

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



for key, grid in vdw_grid.atomic_grids.items():
    charges = fftb.LeastSquaresCharges(molecule, grid)
    charges.sites_to_solution()
    report = fftb.Report(key, charges)
    print report
print molecule
print molecule.frames[0].center.multipoles

