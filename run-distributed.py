import fftoolbox as fftb

import logging, sys

logger = logging.getLogger(__name__)
lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
logging.basicConfig(filename='./log', filemode='w', level=logging.DEBUG, format=lf)
#logging.basicConfig(level=logging.DEBUG, format=lf)

d = {
        'O': (2, 1.),
        'N': (2, 1.),
        'S': (2, 1.),
        'C': (1, 1.),
        'H': (1, 1.),
        }
        

data = {
        'name': 'cis-mesno',
        'theory': 'b3lyp_augccpvdz',
        'representation': ('spherical', 2),
        'symmetry': False,
        'sphere params': d,
        'exclude': ['<xy'],
        }

parser = fftb.GDMA(data=data)
data.update(parser.data)
molecule = fftb.DistributedLebedevMolecule(data=data)


data['density'] = 1.5
parser = fftb.GaussianCube(data=data)
data.update(parser.data.copy())

grid = fftb.vdwGrid(data)
#grid.build_LEM('full-vdw.pymol')
charges = fftb.LeastSquaresCharges(molecule, grid)
charges.sites_to_solution()
report = fftb.Report('full', charges)
print report

atoms = [('O', [1]), 
         ('N', [2]),
         ('S', [3]),
         ('CH3', [4,5,6,7])]

for s, i in atoms:
    data['vdw atoms'] = i
    grid = fftb.vdwGrid(data)
    #grid.build_LEM('atom-%s.pymol' % s)
    charges = fftb.LeastSquaresCharges(molecule, grid)
    charges.sites_to_solution()
    report = fftb.Report(s, charges)
    print report


