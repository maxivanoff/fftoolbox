import numpy as np
import logging, sys
import numpy
import units

__author__ = "Maxim Ivanov"
__email__ = "maxim.ivanov@marquette.edu"

logger = logging.getLogger(__name__)

class LeastSquaresBasic(object):
    """
    LS problem: min ||b-Ax||
    """
    def __init__(self, name=None, A=None, b=None):
        self.name = name
        self.A = A
        self.b = b
        self.At = np.transpose(A)
        self.Atb = np.dot(self.At, self.b)
        self.AtA = np.dot(self.At, self.A)
        self.AAt = np.dot(self.A, self.At)
        self._solution = None
        logger.debug('Shape of reference vector b: %s\nShape of LS matrix A: %s' % (self.b.shape, self.A.shape))

    @property
    def solution(self):
        return self._solution
    
    @property
    def solution_rmsd(self):
        if self.solution is None:
            raise ValueError("self.solution should not be None")
        return float(np.sqrt(self.cost_function(self.solution)))*units.au_to_kcal

    def cost_function(self, x):
        residual = self.b - np.dot(self.A, x)
        return np.average(residual**2)

    def solve(self, method='svd', truncate=0):
        if method == 'normal':
            self._solution = self.solve_normal()
        if method == 'svd':
            self._solution = self.solve_svd(truncate)
        return self._solution

    def solve_normal(self):
        return np.linalg.solve(self.AtA, self.Atb)
    
    def solve_svd(self, truncate_by=0):
        dim = self.A.shape[1]
        self.svdV, self.S, self.svdU = np.linalg.svd(self.A) 
        bv = np.dot(self.b, self.svdV)
        c = bv[:dim]/self.S[:dim]
        c[dim - truncate_by:dim] = np.zeros(truncate_by)
        solution = np.dot(c, self.svdU)
        return solution



