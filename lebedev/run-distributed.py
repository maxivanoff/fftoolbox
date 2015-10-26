from fftoolbox.parser import GaussianCube, QChem, Gaussian, ForceFieldXML, GDMA
from fftoolbox.grid import vdwGrid
from fftoolbox.charges import LeastSquaresCharges
from molecule import DistributedLebedevMolecule
from fftoolbox.results import Results
from charges import LebedevCharges

import logging, sys

#logger = logging.getLogger(__name__)
#lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
#logging.basicConfig(level=logging.DEBUG, format=lf)

d = {
        'H': (3, 1.),
        'N': (3, 1.),
        'O': (3, 1.),
        }

d = {
        'O': (4, 1.),
        'N': (4, 1.),
        'S': (4, 1.),
        'C': (1, 1.),
        'H': (1, 1.),
        }
        

data = {
        'name': 'cis-mesno',
        'theory': 'b3lyp_augccpvdz',
        'representation': ('spherical', 2),
        'symmetry': False,
        'distributions': d,
        }

parser = GDMA(data=data)
#results = Results(parser.data['multipoles'])

data.update(parser.data)

molecule = DistributedLebedevMolecule(data=data)

molecule.write_xyz('mesno-distr.xyz', here=True)

for a in molecule.atoms:
    multipoles = a.multipoles.copy()
    for mname in sorted(multipoles.keys()):
        try:
            a.ref_multipoles[mname]
            print a.name, mname, multipoles[mname], a.ref_multipoles[mname]
        except KeyError:
            print a.name, mname, multipoles[mname], 0.0

multipoles = molecule.multipoles.copy()
for mname in sorted(multipoles.keys()):
    try:
        parser.data['multipoles'][mname]
    except KeyError:
        parser.data['multipoles'][mname] = 0.0
    print mname, multipoles[mname], parser.data['multipoles'][mname]

data['density'] = 1.5
parser = GaussianCube(data=data)
data.update(parser.data.copy())

grid = vdwGrid(data)
grid.build_LEM('full-vdw.pymol')
solution = LebedevCharges(molecule, grid)
print 'full', solution.solution_rmsd

atoms = [('O', [1]), 
         ('N', [2]),
         ('S', [3]),
         ('CH3', [4,5,6,7])]
for s, i in atoms:
    data['vdw atoms'] = i
    grid = vdwGrid(data)
    grid.build_LEM('atom-%s.pymol' % s)
    solution = LebedevCharges(molecule, grid)
    print s, solution.solution_rmsd


