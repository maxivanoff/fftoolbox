import os
import numpy as np
from lebedev_write import Lebedev

WORKDIR = os.path.dirname(os.path.dirname(__file__))
class Sphere(object):
    
    def __init__(self, data=None):
        filename = '%s/data/sphere/%s_%s_p%ir%.1f.qesp' % (WORKDIR, data['name'], data['theory'], data['grid points'], data['grid radius'])
        self.read_file(filename)

    def read_file(self, filename):
        file = open(filename, 'r')
        self.xyz = []
        self.values = []
        self.atoms = []
        self.property = 'esp'
        for line in file:
            tmp = [float(t) for t in line.split()]
            xyz = np.array(tmp[0:3])
            self.xyz.append(xyz)
            self.values.append(tmp[3])
        file.close()

        num_points = len(self.values)
        self.weights = np.array([np.sqrt(4*np.pi*w) for x, y, z, w in Lebedev(num_points)])

        self.data = {
                    'values': self.values,
                    'grid points': self.xyz,
                    'weights': self.weights,
                }

    def __str__(self):
        return "Sphere parser"


