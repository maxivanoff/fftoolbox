import logging, sys
import numpy as np
import cma

from parser import GaussianCube, QChem, Gaussian
from grid import Grid
from units import au_to_kcal
from molecule import Molecule
from charges import LeastSquaresCharges

__author__ = "Maxim Ivanov"
__email__ = "maxim.ivanov@marquette.edu"

class NLLSCharges(LeastSquaresCharges):
    """
    Non Linear Least Squares Charges: positions of extra point charges are optimized
    """

    def __init__(self, molecule=None, grid=None, truncate=0):
        LeastSquaresCharges.__init__(self, molecule=molecule, grid=grid)
        self.truncate = truncate
        self.hybrids = []
        for atom in self.molecule.atoms:
            if not atom.frame is None:
                self.hybrids.append(atom)
        self.num_charges = len(self.molecule.sites_noneq)

    def nonlinear_cost(self, x):
        q = x[0:self.num_charges]
        geom = x[self.num_charges:]
        self.update_geometry(geom)
        return self.cost_function(q)

    def rmsd(self, x):
        return np.sqrt(self.nonlinear_cost(x))*au_to_kcal

    def update_geometry(self, geom):
        for i, atom in enumerate(self.hybrids):
            d, a = geom[i*2:i*2 + 2]
            atom.set_hybridization((atom.frame.hybridization, d, a))
        self.setA()

    def update_charges(self, x):
        charges = {}
        for i, site1 in enumerate(self.molecule.sites_noneq):
            for site2 in self.molecule.get_sites_by_name(site1.name):
                site2.charge = x[i]

    def geom_rmsd(self, geom):
        self.update_geometry(geom)
        self.solve(truncate=self.truncate)
        return self.solution_rmsd

