import fftoolbox as fftb

import logging, sys

logger = logging.getLogger(__name__)
lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
logging.basicConfig(filename='./log', filemode='w', level=logging.DEBUG, format=lf)
#logging.basicConfig(level=logging.DEBUG, format=lf)

d = {
        'C': (0, 0.5),
        'F': (0, 0.5),
        'S': (2, 0.5),
        'H': (0, 0.5),
        }

data = {
        'name': 'cf3sh',
        'theory': 'mp2_augccpvtz',
        'representation': ('spherical', 2),
        'sphere params': d,
        }

parser = fftb.GDMA(data=data)
data.update(parser.data)
molecule = fftb.DistributedLebedevMolecule(data=data)
results = fftb.Results(data['multipoles'])

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

for key, grid in vdw_grid.atomic_grids.items():
    charges = fftb.LeastSquaresCharges(molecule, grid)
    charges.sites_to_solution()
    report = fftb.Report(key, charges)
    print report

