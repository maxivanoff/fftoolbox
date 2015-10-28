from fftoolbox.charges import LeastSquaresCharges
import numpy as np

class LebedevCharges(LeastSquaresCharges):

    def __init__(self, molecule=None, grid=None):
        LeastSquaresCharges.__init__(self, molecule=molecule, grid=grid)

    @property
    def solution(self):
        self._solution = np.zeros(self.molecule.num_sites)
        for i, s in enumerate(self.molecule.sites):
            self._solution[i] = s.charge
        return self._solution
        

