from fftoolbox2.parser import GaussianCube, QChem, Gaussian
from fftoolbox2.grid import vdwGrid
from fftoolbox2.constrained import ConstrainedLeastSquaresCharges, LagrangeConstraint, TrivialConstraint, SVDConstraint
from fftoolbox2.charges import LeastSquaresCharges
from fftoolbox2.molecule import Molecule
from fftoolbox2.results import Results
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

import logging, sys

logger = logging.getLogger(__name__)
lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
logging.basicConfig(level=logging.DEBUG, format=lf)

data = {
        'name': 'methanol',
        'theory': 'b3lyp_augccpvdz',
        'density': 1.5,
        'exclude': [],
        'include': [],
        'symmetry': False,
        'multipoles': ('cartesian', 2)
        }

parser = GaussianCube(data=data)
#results = Results(Gaussian(data=data).data['multipoles'])
results = Results()

data.update(parser.data.copy())

grid = vdwGrid(data)
grid.build_LEM(filename='./%s.pymol' % data['name'])
molecule = Molecule(data)


ls = LeastSquaresCharges(grid=grid, molecule=molecule)
ls.solve()
results.add(ls, name='no contraint')

ls = LagrangeConstraint(grid=grid, molecule=molecule)
ls.solve()
results.add(ls, name='lagrange')

ls = TrivialConstraint(grid=grid, molecule=molecule)
ls.solve()
results.add(ls, name='trivial')

ls = SVDConstraint(grid=grid, molecule=molecule)
ls.solve(truncate=0)
results.add(ls, name='svd 0')

ls = SVDConstraint(grid=grid, molecule=molecule)
ls.solve(truncate=1)
results.add(ls, name='svd 1')

for sname in molecule.sites_names_noneq:
    ls = ConstrainedLeastSquaresCharges(grid=grid, molecule=molecule, eliminated=sname) 
    ls.solve(method='svd', truncate=0)
    results.add(ls, name=sname)

print results

