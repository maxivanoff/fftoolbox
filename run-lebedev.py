import fftoolbox as fftb

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
parser = fftb.GaussianCube(data=data)
data.update(parser.data.copy())
grid = fftb.vdwGrid(data)

for rank in xrange(7):
    data = {
            'name': 'cis-mesno',
            'theory': 'b3lyp_augccpvdz',
            'sphere params': (rank, 1.0),
            }

    parser = fftb.GDMA(data=data)
    data.update(parser.data.copy())
    molecule = fftb.LebedevMolecule(data)

    #grid.build_LEM('full-vdw.pymol')
    charges = fftb.LeastSquaresCharges(molecule, grid)
    charges.sites_to_solution()
    report = fftb.Report('rank %i' % rank, charges)
    print report
 
