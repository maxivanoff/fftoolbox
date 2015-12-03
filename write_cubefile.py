import fftoolbox as fftb

import logging, sys

logger = logging.getLogger(__name__)
lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
logging.basicConfig(level=logging.DEBUG, format=lf)

data = {
        'name': 'methanethiol',
        'theory': 'mp2_augccpvtz',
        }

parser = fftb.GaussianCube(data=data)
data.update(parser.data.copy())

cubic_grid = fftb.Grid(data)
molecule = fftb.HM(data)

writer = fftb.GaussianCube()
qm_values = parser.data['values']*fftb.au_to_kcal
writer.write_file(grid=cubic_grid, molecule=molecule,\
                  values=qm_values, filename='%s.cub' % (data['name']))

data = {
        'name': 'methanethiol',
        'theory': 'mp2_augccpvtz',
        'density': 1.5,
        'symmetry': False,
        'exclude': ['<xy'],
        'representation': ('cartesian', 1)
        }

parser = fftb.GaussianCube(data=data)
data.update(parser.data.copy())
vdw_grid = fftb.vdwGrid(data)

ls = fftb.LSC(grid=vdw_grid, molecule=molecule) 
ls.solve()

cubic_grid.create_grid()
mep = fftb.MEP(cubic_grid, molecule)
ac_values = mep.compute('point charges')

writer.write_file(grid=cubic_grid, molecule=molecule,\
                  values=ac_values, filename='%s_ac.cub' % (data['name']))

