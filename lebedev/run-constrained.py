from parser import Sphere
from grid import LebedevGrid
from constrained import LebedevConstrained
from charges import LebedevLeastSquares
from molecule import LebedevMolecule
from fftoolbox2.results import Results
import fftoolbox2.visual as vis
import numpy as np

import logging, sys

logger = logging.getLogger(__name__)
lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
logging.basicConfig(level=logging.DEBUG, format=lf)

data = {
        'name': 'methanol',
        'theory': 'b3lyp_augccpvdz',
        'symmetry': False,
        'multipoles': ('spherical', 3),
        'grid points': 194,
        'grid radius': 8.0,
        'points': 6,
        'radius': 2.0,
        }

results = Results()
parser = Sphere(data=data)
data.update(parser.data.copy())

grid = LebedevGrid(data)
grid.build_LEM(filename='./lebedev.pymol' )
molecule = LebedevMolecule(data)

ls = LebedevLeastSquares(grid=grid, molecule=molecule) 
ls.solve()
results.add(ls)

U = np.transpose(ls.svdU)
vis.table(U, 'methanol', molecule.sites_names_noneq, ls.S)

ls = LebedevConstrained(grid=grid, molecule=molecule, eliminated='2') 
ls.solve()
results.add(ls)

sites_names = [s for s in molecule.sites_names_noneq if not s=='2']
U = np.transpose(ls.svdU)
vis.table(U, '2', sites_names, ls.S)

print results

