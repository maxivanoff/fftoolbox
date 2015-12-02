import numpy as np
import logging

from LeastSquares import LeastSquaresBasic
from molecule import Molecule
from grid import Grid
from atom import Coordinates
import fast
import units

__author__ = "Maxim Ivanov"
__email__ = "maxim.ivanov@marquette.edu"


logger = logging.getLogger(__name__)

class LeastSquaresCharges(LeastSquaresBasic):

    def __init__(self, molecule=None, grid=None):
        self.molecule = molecule
        self.update_grid(grid)
        self._charges = {}

    def update_grid(self, grid):
        self.grid = grid
        self.setA()
        LeastSquaresBasic.__init__(self, A=self.A, b=grid.get_values())

    def setA_slow(self):
        n_p = len(self.grid.points)
        n_s = len(self.molecule.sites_noneq)
        self.A =  np.zeros((n_p, n_s))
        for i in xrange(n_p):
            proton = Coordinates(coordinates=self.grid.points[i].coordinates)
            for j, name in enumerate(self.molecule.sites_names_noneq):
                for site in self.molecule.get_sites(name=name):
                    self.A[i,j] += 1/site.distance_to(proton)

    def setA(self):
        grid_coordinates = self.grid.get_coordinates()
        sites_coordinates = self.molecule.get_coordinates()
        num_reduced_sites = len(self.molecule.sites_names_noneq) 
        logger.debug('Number of coordinates in the grid: %i\nNumber of sites in molecule: %i\nNumber of sites with different names: %i' % (len(grid_coordinates), len(sites_coordinates), num_reduced_sites))
        self.A = fast.set_inversed(grid_coordinates, sites_coordinates, \
                num_reduced_sites, self.molecule.sym_sites)
        logger.info("Least squares matrix is set up with shape (%s, %s)" % (self.A.shape[0], self.A.shape[1]))

    @property
    def charges(self):
        return self._charges

    def clear_charges(self):
        self._charges = {}

    def set_charge(self, name, q):
        self._charges[name] = q

    @property
    def multipoles(self):
        return self.molecule.charges_to_multipoles(self.charges)

    def solve(self, method='svd', truncate=0):
        solution = LeastSquaresBasic.solve(self, method, truncate)
        if not len(solution) == len(self.molecule.sites_names_noneq):
            raise ValueError('Number of charges in solution (%i) does not match number of sites (%i):\n%r' % 
                    (len(solution), len(self.molecule.sites_names_noneq), self.molecule.sites_names_noneq))
        self.clear_charges()
        for q, name in zip(solution, self.molecule.sites_names_noneq):
            self.set_charge(name, q)
        self.charges_to_sites(self.charges)

    def sites_to_solution(self):
        solution = np.zeros(self.molecule.num_sites_noneq)
        for i, s in enumerate(self.molecule.sites_noneq):
            solution[i] = s.charge
        self.set_solution(solution)
        logger.info('Charges were transfered from sites to _solution')

    def charges_to_sites(self, charges):
        for site in self.molecule.sites:
            q = charges[site.name]
            site.set_charge(q)
        logger.info('Charges were transfered to sites')


