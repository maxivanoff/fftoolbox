from scipy.special import sph_harm as Y
from fftoolbox.multipole import Ylmc, Ylms
from numpy.linalg import norm
import numpy as np
import units
from time import time

class MEP(object):

    Z = {'H': 1, 'C': 4, 'N': 5, 'O':6, 'S':16}

    def __init__(self, grid=None, molecule=None):
        self.grid = grid
        self.molecule = molecule

    def compute(self, representation='point charges'):
        if representation == 'point charges':
            values = self.point_charge()
        if representation == 'multipole':
            values = self.multipole()
        return values

    def point_charge(self):
        values = np.zeros(self.grid.num_points)
        for i, p in enumerate(self.grid.points):
            f = 0.
            for s in self.molecule.sites:
                r = np.linalg.norm(p.coordinates - s.coordinates)
                if abs(r)< 0.0001:
                    f += 0.
                else:
                    f += s.charge/r
            values[i] = f
        values *= units.au_to_kcal
        return values

    def multipole(self):
        values = np.zeros(self.grid.num_points)
        for i, p in enumerate(self.grid.points):
            f = 0.
            for name, mult_data in self.molecule.multipoles_data().items():
                origin, rank, multipoles = mult_data
                r, theta, phi = p.cartesian_to_spherical(origin)
                for l in range(rank+1):
                    if abs(r) < 0.0001:
                        k = 0.
                    else:
                        k = np.sqrt(4*np.pi/(2*l + 1))*np.power(r, -l-1)
                    for m in range(-l, l+1):
                        if m < 0:
                            mname = '%i%is' % (l, abs(m))
                            YY = Ylms
                        if m == 0:
                            mname = '%i0' % l
                            YY = Y
                        if m > 0:
                            mname = '%i%ic' % (l, m)
                            YY = Ylmc
                        try:
                            Q = multipoles[mname]
                        except KeyError:
                            Q = 0.
                        f += k*Q*YY(abs(m), l, theta, phi).real
            values[i] = f
        values *= units.au_to_kcal
        return values



