import numpy as np
from fftoolbox.charges import LeastSquaresCharges
import logging
logger = logging.getLogger(__name__)

class LebedevLeastSquares(LeastSquaresCharges):

    def __init__(self, molecule=None, grid=None):
        LeastSquaresCharges.__init__(self, molecule=molecule, grid=grid)
        self.updateA()

    def updateA(self):
        for i, wp in enumerate(self.grid.weights):
            for j, ws in enumerate(self.molecule.weights):
                self.A[i,j] *= wp*ws

    def cost_function(self, x):
        residual = self.b - np.dot(self.A, x)
        residual /= self.grid.weights
        return np.average(residual**2)

    @property
    def charges(self):
        logger.info('Charges were multiplied by quadrature weights')
        charges = {}
        for sigma, name, w in zip(self.solution, self.molecule.sites_names, self.molecule.weights):
            charges[name] = sigma*w
        return charges



