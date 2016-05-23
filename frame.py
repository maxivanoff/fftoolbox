from __future__ import division
from numpy.linalg import norm
import numpy as np
import logging
import units

__author__ = "Maxim Ivanov"
__email__ = "maxim.ivanov@marquette.edu"

logger = logging.getLogger(__name__)

class BasicFrame(object):

    def __init__(self, frame=None):
        self.local_axes = np.identity(3)
        if frame:
            a, b, c = frame
            self.set_axes(a, b, c)

    def set_axes(self, a, b, c):
        ab = a - b
        ac = a - c
        ab_normed, ac_normed = [v/norm(v) for v in [ab, ac]]
        x = ab_normed + ac_normed
        y = ab_normed - ac_normed
        z = np.cross(x, y)
        for i, a in enumerate([x,y,z]):
            self.local_axes[:,i] = a/norm(a)
        logger.debug("Local coordinate system is set:\n%r" % self.local_axes)

class AtomFrame(BasicFrame):

    def __init__(self, atom):
        self.center = atom
        self.neighbors = atom.neighbors
        self.num_nghbrs = len(self.center.neighbors)
        logger.info('Frame center at %s has %i neighbors' % (self.center.name, self.num_nghbrs))
        frame = self.get_frame_vectors()
        BasicFrame.__init__(self, frame)

    def get_frame_vectors(self):
        a = self.center.coordinates
        if self.num_nghbrs == 2:# e.g. OH2 or HN=CH2
            b = self.center.neighbors[0].coordinates
            c = self.center.neighbors[1].coordinates
            logger.info('a is %s, b is %s, c is %s' % (self.center.name, self.center.neighbors[0].name, self.center.neighbors[1].name))
        if self.num_nghbrs == 1: # e.g. O=CH2
            X = self.center.neighbors[0]
            if len(X.neighbors) == 3:
                x1, x2 = [atom.coordinates for atom in X.neighbors if not atom.name == self.center.name]
                b = (x1 + X.coordinates)*0.5
                c = (x2 + X.coordinates)*0.5
            if len(X.neighbors) == 2: # e.g. O=N-S-
                x1 = [atom.coordinates for atom in X.neighbors if not atom.name == self.center.name][0]
                c = (x1 - X.coordinates)/norm(x1 - X.coordinates) + X.coordinates
                v1 = X.coordinates - X.neighbors[0].coordinates
                v2 = X.coordinates - X.neighbors[1].coordinates
                b = X.coordinates + (v1 + v2)/norm(v1 + v2)
        return a, b, c

