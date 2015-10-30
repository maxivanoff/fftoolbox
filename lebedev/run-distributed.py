from fftoolbox.parser import GaussianCube, QChem, Gaussian, ForceFieldXML, GDMA
from fftoolbox.grid import vdwGrid
from fftoolbox.charges import LeastSquaresCharges
from molecule import DistributedLebedevMolecule
from fftoolbox.results import Results
from charges import LebedevCharges

import logging, sys

logger = logging.getLogger(__name__)
lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
logging.basicConfig(filename='./log', filemode='w', level=logging.DEBUG, format=lf)
#logging.basicConfig(level=logging.DEBUG, format=lf)

d = {
        'H': (3, 1.),
        'N': (3, 1.),
        'O': (3, 1.),
        }

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

parser = GDMA(data=data)
#results = Results(parser.data['multipoles'])
data.update(parser.data)
molecule = DistributedLebedevMolecule(data=data)


data['density'] = 1.5
parser = GaussianCube(data=data)
data.update(parser.data.copy())

grid = vdwGrid(data)
#grid.build_LEM('full-vdw.pymol')
charges = LeastSquaresCharges(molecule, grid)
charges.sites_to_solution()
print 'full', charges.rmsd

atoms = [('O', [1]), 
         ('N', [2]),
         ('S', [3]),
         ('CH3', [4,5,6,7])]
for s, i in atoms:
    data['vdw atoms'] = i
    grid = vdwGrid(data)
    #grid.build_LEM('atom-%s.pymol' % s)
    charges = LeastSquaresCharges(molecule, grid)
    charges.sites_to_solution()
    print s, charges.rmsd


