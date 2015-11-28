import fftoolbox as fftb

import logging, sys

logger = logging.getLogger(__name__)
lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
#logging.basicConfig(level=logging.DEBUG, format=lf)
logging.basicConfig(filename='./log', filemode='w', level=logging.DEBUG, format=lf)

data = {
        'name': 'methanethiol',
        'theory': 'mp2_augccpvtz',
        'density': 1.5,
        }
parser = fftb.GaussianCube(data=data)
data.update(parser.data.copy())
grid = fftb.vdwGrid(data)

radius = 0.5
for rank in xrange(1,4):
    data = {
            'name': 'methanethiol',
            'theory': 'mp2_augccpvtz',
            'sphere params': (rank, radius),
            }

    parser = fftb.GDMA(data=data)
    data.update(parser.data.copy())
    molecule = fftb.LebedevMolecule(data)
    molecule.color_charges(filename='%s_%s-%i-%.1f.pymol' % (data['name'], data['theory'], rank, radius), xyzname='%s_%s.xyz' % (data['name'], data['theory']))

    #grid.build_LEM('full-vdw.pymol')
    charges = fftb.LeastSquaresCharges(molecule, grid)
    charges.sites_to_solution()
    report = fftb.Report('rank %i' % rank, charges)
    print report
 
