import numpy as np
import logging

from LeastSquares import LeastSquaresBasic
from molecule import Molecule
from grid import Grid
from atom import Site
import fast
import units

__author__ = "Maxim Ivanov"
__email__ = "maxim.ivanov@marquette.edu"


logger = logging.getLogger(__name__)

class LeastSquaresCharges(LeastSquaresBasic):

    def __init__(self, molecule=None, grid=None):
        self.molecule = molecule
        self.grid = grid
        self.setA()
        LeastSquaresBasic.__init__(self, A=self.A, b=grid.get_values())

    def setA_slow(self):
        n_p = len(self.grid.points)
        n_s = len(self.molecule.sites_noneq)
        self.A =  np.zeros((n_p, n_s))
        for i in xrange(n_p):
            proton = Site(coordinates=self.grid.points[i].coordinates, name='H+')
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
        charges = {}
        for q, name in zip(self.solution, self.molecule.sites_names):
            charges[name] = q
        return charges

    @property
    def multipoles(self):
        return self.molecule.charges_to_multipoles(self.charges)

    def charges2sites(self):
        current_charges = self.charges.copy()
        for i, site in enumerate(self.molecule.sites):
            site.charge = current_charges[site.name]
        logger.info('Charges were transfered to sites')


