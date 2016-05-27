import fftoolbox as fftb

import logging, sys

#logger = logging.getLogger(__name__)
#lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
#logging.basicConfig(level=logging.DEBUG, format=lf)
#logging.basicConfig(filename='./log', filemode='w', level=logging.DEBUG, format=lf)

molecule_name = 'sf2'
radius = 0.5
rank = 3

data = {
        'name':molecule_name,
        'theory': 'b3lyp_augccpvdz',
        'density': 1.5,
        }
parser = fftb.GaussianCube(data=data)
data.update(parser.data.copy())
grid = fftb.vdwGrid(data)

data = {
        'name': molecule_name,
        'theory': 'b3lyp_augccpvdz',
        'sphere params': (rank, radius),
        }

parser = fftb.GDMA(data=data)
data.update(parser.data.copy())
molecule = fftb.LebedevMolecule(data)

molecule.write_mdgx_in('SF2', 8, 'S', 'F1', 'F2', fname='%s-l%i.xpt' % (molecule_name, rank))
    
charges = fftb.LeastSquaresCharges(molecule, grid)
for s in molecule.sites:
    charges.set_charge(s.name, s.charge)
results = fftb.Results(data['multipoles'])
results.add(charges)
print results

 
