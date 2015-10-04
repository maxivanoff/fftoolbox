from parser import Sphere
from grid import LebedevGrid
from charges import LebedevLeastSquares
from molecule import LebedevMolecule
from fftoolbox.results import Results

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

print results

