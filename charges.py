import units
from LeastSquares import LeastSquaresBasic
from molecule import Molecule
from grid import Grid
from atom import Site
import numpy as np
import fast
import time

import logging

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
        self.A = fast.set_inversed(grid_coordinates, sites_coordinates, \
                num_reduced_sites, self.molecule.sym_sites)

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


