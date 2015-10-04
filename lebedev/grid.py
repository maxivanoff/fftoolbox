from fftoolbox.grid import GridPoint, Grid
import numpy as np
import logging
import os

logger = logging.getLogger(__name__)

class LebedevGrid(Grid):

    def __init__(self, data):
        self.molecule_name = data['name']
        self.theory = data['theory']
        self.points = [] 
        self.radius = data['grid radius']
        self.create_grid(points=data['grid points'], values=data['values'])
        self.weights = data['weights']

    def create_grid(self, points, values):
        self.points = []
        logger.info("Creating grid from list of values.\nNumber of grid points: %i\nNumber of grid values: %i\nRadius: %.1f" % (len(points), len(values), self.radius))
        for p, v in zip(points, values):
            point = GridPoint(coordinates=p, value=v)
            self.points.append(point)

    def get_values(self):
        logger.info("Grid point values were multiplied by the quadrature weights\nNumber of weights:%i\nNumber of sites:%i" % (len(self.weights), len(self.points)))
        return np.array([p.value*w for p, w in zip(self.points, self.weights)])
