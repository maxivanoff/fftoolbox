from scipy.optimize import curve_fit
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
        self.set_b(b)
        self.set_A(A)
        self._solution = None
        logger.debug('Shape of reference vector b: %s\nShape of LS matrix A: %s' % (self.b.shape, self.A.shape))

    def set_A(self, A):
        self.A = A
        self.At = np.transpose(A)
        self.Atb = np.dot(self.At, self.b)
        self.AtA = np.dot(self.At, self.A)
        self.AAt = np.dot(self.A, self.At)

    def set_b(self, b):
        self.b = b

    def set_solution(self, solution):
        self._solution = solution

    @property
    def solution(self):
        return self._solution
    
    def copy_solution(self):
        return self._solution.copy()

    @property
    def residual(self):
        return self.b - np.dot(self.A, self.solution)
    
    @property
    def rmsd(self):
        rmsd2 = np.average(self.residual**2)
        return np.sqrt(rmsd2)*units.au_to_kcal

    @property
    def R2(self):
        y = np.dot(self.A, self.solution)
        correlation = np.corrcoef(self.b, y)[0,1]
        return correlation**2

    @property
    def max_error(self):
        return np.amax(np.absolute(self.residual))*units.au_to_kcal

    @property
    def rmad(self):
        return np.sum(np.absolute(self.residual))/np.sum(np.absolute(self.b))

    @property
    def ab(self):
        x = self.b
        y = np.dot(self.A, self.solution)
        popt, pcov = curve_fit(lambda xdata,m,n: m*xdata+n, x, y)
        a, b = popt
        da, db = np.sqrt(np.diag(pcov))
        return a, da, b, db

    def cost_function(self, x):
        residual = self.b - np.dot(self.A, x)
        return np.average(residual**2)

    def solve(self, method='svd', truncate=0):
        if method == 'normal':
            logger.info('Solving least squares problem using normal equations')
            s = self.solve_normal()
        if method == 'svd':
            logger.info('Solving least squares problem using SVD and truncating %i terms' % truncate)
            s = self.solve_svd(truncate)
        self.set_solution(s)
        return self.solution

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



