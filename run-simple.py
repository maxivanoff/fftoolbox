from parser import GaussianCube, QChem, Gaussian, ForceFieldXML
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
        'representation': ('cartesian', 2)
        }

results = Results(Gaussian(data=data).data['multipoles'])

parser = GaussianCube(data=data)
data.update(parser.data.copy())

grid = vdwGrid(data)
molecule = Molecule(data)

ls = LeastSquaresCharges(grid=grid, molecule=molecule) 
ls.solve()
ls.charges2sites()
results.add(ls)

print results

# save
molecule.write_xyz(filename='%s_%s.xyz' % (data['name'], data['name']))
molecule.write_mol2(filename='%s_%s.mol2' % (data['name'], data['name']))

ff = ForceFieldXML()
ff.write_file(molecule=molecule)
ff.load_forcefields(filename='%s_%s.xml' % (data['name'], data['theory']), molecule=molecule)


