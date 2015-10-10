from parser import GaussianCube, QChem, Gaussian, ForceFieldXML
from grid import vdwGrid
from charges import LeastSquaresCharges
from molecule import Molecule
from results import Results

import logging, sys

logger = logging.getLogger(__name__)
lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
logging.basicConfig(level=logging.DEBUG, format=lf)

amber = {
        'O-sp3':('sp3', 0.35, 109.5/2.),
        'S-sp3':('sp3', 0.70, 90.),
        'O-sp2': ('sp2', 0.35, 60.0),
        'N-sp2': ('sp2', 0.35, 0.0)
        }

data = {
        'name': 'methanol',
        'theory': 'pbe_def2svp',
        'density': 1.5,
        'symmetry': False,
        'hybridizations': {'O': amber['O-sp3']},
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
molecule.write_xyz(filename='amber-%s_%s.xyz' % (data['name'], data['theory']))
molecule.write_mol2(filename='amber-%s_%s.mol2' % (data['name'], data['theory']))

ff = ForceFieldXML()
ff.write_file(molecule=molecule, xmlfilename='data/forcefields/amber-%s_%s.xml' % (data['name'], data['theory']))
ff.load_forcefields(filename='amber-%s_%s.xml' % (data['name'], data['theory']), molecule=molecule)


