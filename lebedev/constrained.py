from fftoolbox.atom import Site
from fftoolbox.constrained import ConstrainedLeastSquaresCharges
import numpy as np

import logging

logger = logging.getLogger(__name__)

class LebedevConstrained(ConstrainedLeastSquaresCharges):

    def __init__(self, molecule=None, grid=None, eliminated=None):
        ConstrainedLeastSquaresCharges.__init__(self, molecule=molecule, grid=grid, eliminated=eliminated)
        self.updateA()

    def updateA(self):
        reduced_weights = [w for w, s in zip(self.molecule.weights, self.molecule.sites) if not s.name == self.eliminated]
        for i, wp in enumerate(self.grid.weights):
            for j, ws in enumerate(reduced_weights):
                self.A[i,j] *= wp*ws

    def cost_function(self, x):
        residual = self.b - np.dot(self.A, x)
        residual /= self.grid.weights
        return np.average(residual**2)

    @property
    def charges(self):
        charges = {}
        sites = [s.name for s in self.molecule.sites if not s.name == self.eliminated]
        reduced_weights = [w for w, s in zip(self.molecule.weights, self.molecule.sites) if not s.name == self.eliminated]
        solution = []
        for sigma, name, w in zip(self.solution, sites, reduced_weights):
            charges[name] = sigma*w
            solution.append(sigma*w)
        charges[self.eliminated] = -np.sum(solution)
        return charges



