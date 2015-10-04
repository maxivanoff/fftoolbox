from numpy.linalg import norm
import numpy as np
import units

class MEP(object):

    Z = {'H': 1, 'C': 4, 'N': 5, 'O':6, 'S':16}

    def __init__(self, grid=None, molecule=None):
        self.grid = grid
        self.molecule = molecule
        self.evaluate()

    def evaluate(self):
        values = []
        for p in self.grid.points:
            f = 0.
            for s in self.molecule.sites:
                r = norm(p.coordinates - s.coordinates)
                #f += (1 - np.exp(-2*r)*(1+r))*s.charge/r
                f += s.charge/r
                if not s.element[0] == 'E':
                    Z = self.Z[s.element]
                    f += np.exp(-2.*r)*(Z - (s.charge - Z)/r)
            values.append(f)
        self.values = np.array(values)*units.au_to_kcal

