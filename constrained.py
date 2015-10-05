import numpy as np
import logging

from atom import Site
from charges import LeastSquaresCharges
import visual

__author__ = "Maxim Ivanov"
__email__ = "maxim.ivanov@marquette.edu"

logger = logging.getLogger(__name__)

class ConstrainedLeastSquaresCharges(LeastSquaresCharges):

    def __init__(self, molecule=None, grid=None, eliminated=None):
        self.eliminated = eliminated
        logger.info('Applying total charge constraint by elimination.\nEliminated atom: %s' % eliminated)
        LeastSquaresCharges.__init__(self, molecule=molecule, grid=grid)
        self.setB()
        self.A -= self.B

    def setA(self):
        n_p = len(self.grid.points)
        n_s = len(self.molecule.sites) - 1
        self.A =  np.zeros((n_p, n_s))
        for i in xrange(n_p):
            proton = Site(coordinates=self.grid.points[i].coordinates, name='H+')
            for j, site in enumerate([s for s in self.molecule.sites if not s.name == self.eliminated]):
                self.A[i,j] += 1/site.distance_to(proton)

    def setB(self):
        n_p = len(self.grid.points)
        n_s = len(self.molecule.sites) - 1
        self.B =  np.zeros((n_p, n_s))
        X = self.molecule.get_sites(name=self.eliminated)[0]
        for i in xrange(n_p):
            proton = Site(coordinates=self.grid.points[i].coordinates, name='H+')
            rin = 1/X.distance_to(proton)
            for j in xrange(n_s):
                self.B[i,j] = rin

    @property
    def charges(self):
        charges = {}
        sites = [s.name for s in self.molecule.sites if not s.name == self.eliminated]
        for q, name in zip(self.solution, sites):
            charges[name] = q
        charges[self.eliminated] = -np.sum(self.solution)
        return charges

class LagrangeConstraint(LeastSquaresCharges):

    def __init__(self, molecule=None, grid=None):
        LeastSquaresCharges.__init__(self, molecule=molecule, grid=grid)
        logger.info('Applying total charge constraint using Lagrange multiplier')
        self.update()

    def update(self):
        ind = [self.molecule.sites_names_eq.count(name) for name in self.molecule.sites_names_noneq]
        s1, s2 = self.AtA.shape
        AtA = np.zeros((s1+1, s2+1))
        AtA[:-1,:-1] = self.AtA
        AtA[-1,:-1] = ind[:]
        AtA[:-1,-1] = ind[:]
        self.AtA = AtA
        self.Atb = np.append(self.Atb, 0)

    def solve(self):
        solution = LeastSquaresCharges.solve_normal(self)
        self._solution = solution[:-1]
        return self._solution

class TrivialConstraint(LeastSquaresCharges):

    def __init__(self, molecule=None, grid=None):
        LeastSquaresCharges.__init__(self, molecule=molecule, grid=grid)
        logger.info('The total charge leakage will be distributed over all atoms')

    @property
    def solution(self):
        ind = np.array([self.molecule.sites_names_eq.count(name) for name in self.molecule.sites_names_noneq])
        Q = np.sum(ind*self._solution)
        dQ = - np.ones(len(self._solution))*Q/len(self.molecule.sites_names_noneq)
        return self._solution + dQ 

class SVDConstraint(LeastSquaresCharges):

    def __init__(self, molecule=None, grid=None):
        LeastSquaresCharges.__init__(self, molecule=molecule, grid=grid)
        logger.info('Charge constraint using SVD')

    def solve(self, truncate=0):
        v, s, u = np.linalg.svd(self.A)
        u[0,:] = np.ones(u.shape[0])
        c = np.array([np.dot(self.b, v[::,i])/s[i] for i in xrange(u.shape[0])])
        c[0] = 0.
        for i in xrange(truncate):
            c[-1-i] = 0.
        self._solution = np.linalg.solve(u,c)
        logger.debug('Truncated SVD with constraint is used.\nc = %s' % c)
        return self._solution
