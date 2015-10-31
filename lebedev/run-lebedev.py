from fftoolbox.grid import vdwGrid
from fftoolbox.charges import LeastSquaresCharges
from fftoolbox.parser import GDMA, GaussianCube
from molecule import LebedevMolecule
from fftoolbox.results import Report

import logging, sys

logger = logging.getLogger(__name__)
lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
#logging.basicConfig(level=logging.DEBUG, format=lf)
logging.basicConfig(filename='./log', filemode='w', level=logging.DEBUG, format=lf)

data = {
        'name': 'cis-mesno',
        'theory': 'b3lyp_augccpvdz',
        'density': 1.5,
        }
parser = GaussianCube(data=data)
data.update(parser.data.copy())
grid = vdwGrid(data)

for rank in xrange(7):
    data = {
            'name': 'cis-mesno',
            'theory': 'b3lyp_augccpvdz',
            'sphere params': (rank, 1.0),
            }

    parser = GDMA(data=data)
    data.update(parser.data.copy())
    molecule = LebedevMolecule(data)

    #grid.build_LEM('full-vdw.pymol')
    charges = LeastSquaresCharges(molecule, grid)
    charges.sites_to_solution()
    report = Report('rank %i' % rank, charges)
    print report
 
