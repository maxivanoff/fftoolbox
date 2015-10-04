from parser import GaussianCube, QChem, Gaussian
from grid import vdwGrid
from charges import LeastSquaresCharges
from molecule import Molecule
from results import Results

import logging, sys

logger = logging.getLogger(__name__)
lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
logging.basicConfig(level=logging.DEBUG, format=lf)

data = {
        'name': 'methanol',
        'theory': 'pbe_def2svp',
        'density': 1.5,
        'symmetry': False,
        'multipoles': ('cartesian', 2)
        }

results = Results(Gaussian(data=data).data['multipoles'])

parser = GaussianCube(data=data)
data.update(parser.data.copy())

grid = vdwGrid(data)
molecule = Molecule(data)

ls = LeastSquaresCharges(grid=grid, molecule=molecule) 
ls.solve()
results.add(ls)

print results

