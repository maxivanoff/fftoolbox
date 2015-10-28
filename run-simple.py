from parser import GaussianCube, QChem, Gaussian, ForceFieldXML
from grid import vdwGrid
from charges import LeastSquaresCharges
from molecule import HybridMolecule
from results import Results

import logging, sys

logger = logging.getLogger(__name__)
lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
logging.basicConfig(level=logging.DEBUG, format=lf)

data = {
        'name': 'cis-mesno',
        'theory': 'b3lyp_augccpvdz',
        'density': 1.5,
        'symmetry': False,
        'representation': ('cartesian', 1)
        }

results = Results(Gaussian(data=data).data['multipoles'])

parser = GaussianCube(data=data)
data.update(parser.data.copy())

grid = vdwGrid(data)
molecule = HybridMolecule(data)

ls = LeastSquaresCharges(grid=grid, molecule=molecule) 
ls.solve()
results.add(ls)

print results

# save
molecule.write_xyz(filename='%s_%s.xyz' % (data['name'], data['theory']))
molecule.write_mol2(filename='%s_%s.mol2' % (data['name'], data['theory']))

ff = ForceFieldXML()
ff.write_file(molecule=molecule)
ff.load_forcefields(filename='%s_%s.xml' % (data['name'], data['theory']), molecule=molecule)


