import fftoolbox as fftb

import logging, sys

logger = logging.getLogger(__name__)
lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
logging.basicConfig(filename='./log', filemode='w', level=logging.DEBUG, format=lf)
#logging.basicConfig(level=logging.DEBUG, format=lf)

d = {
        'O': (2, 0.5),
        'N': (2, 0.5),
        'S': (2, 0.5),
        'C': (2, 0.5),
        'H': (2, 0.5),
        }
        

data = {
        'name': 'methanethiol',
        'theory': 'mp2_augccpvtz',
        'representation': ('spherical', 2),
        'symmetry': False,
        'sphere params': (2, 0.5),
        'exclude': [],
        }

parser = fftb.GDMA(data=data)
data.update(parser.data)
molecule = fftb.DistributedLebedevMolecule(data=data)
molecule.write_xyz('%s_%s-2-0.5.xyz' % (data['name'], data['theory']), here=True)
molecule.color_charges('%s_%s-2-0.5-DLM.pymol' % (data['name'], data['theory']), xyzname='%s_%s.xyz' % (data['name'], data['theory']))


data['density'] = 1.5
parser = fftb.GaussianCube(data=data)
data.update(parser.data.copy())
grid = fftb.vdwGrid(data)

for key, grid in grid.atomic_grids.items():
    charges = fftb.LeastSquaresCharges(molecule, grid)
    charges.sites_to_solution()
    report = fftb.Report(key, charges)
    print report

