from __future__ import division
from numpy.linalg import norm
import numpy as np
import logging
import units

__author__ = "Maxim Ivanov"
__email__ = "maxim.ivanov@marquette.edu"

logger = logging.getLogger(__name__)

class Frame(object):

    def __init__(self, atom, hybridization):
        self.hybridization = hybridization
        self.center = atom
        self.neighbors = atom.neighbors
        self.num_ep = int(hybridization[-1]) + 1 - len(self.neighbors)
        if hybridization == 'sp3': self.offset = 0. # extra points are in yz plane
        if hybridization == 'sp2': self.offset = np.pi/2. # extra points are in xz plane
        
        logger.info("Creating Frame instance for %s.\nHybridization: %s\nNumber of extra points: %i" % (atom.name, hybridization, self.num_ep))
        self.set_local_system()
        logger.debug("Local coordinate system:\n%r" % self.local_system)
        
    def set_local_system(self):
        a, b, c = self.set_frame()
        ba = b - a
        bc = b - c
        bc_normed, ba_normed = [v/norm(v) for v in [ba, bc]]
        z = bc_normed + ba_normed
        x = bc_normed - ba_normed
        y = np.cross(x,z)
        self.local_system = np.zeros((3,3))
        for i, a in enumerate([x,y,z]):
            self.local_system[:,i] = a/norm(a)

    def set_frame(self):
        b = self.center.coordinates
        if len(self.center.neighbors) == 2: # e.g. OH2 or HN=CH2
            logger.debug('Frame center at %s has 2 neighbors and %i extra points' % (self.center.name, self.num_ep))
            a, c = [atom.coordinates for atom in self.center.neighbors]
        if len(self.center.neighbors) == 1: # e.g. O=CH2
            logger.debug('Frame center at %s has 1 neighbors and %i extra points' % (self.center.name, self.num_ep))
            X = self.center.neighbors[0]
            if len(X.neighbors) == 3:
                x1, x2 = [atom.coordinates for atom in X.neighbors if not atom.name == self.center.name]
                a = (x1 + b)/2.
                c = (x2 + b)/2.
            if len(X.neighbors) == 2:
                x1 = [atom.coordinates for atom in X.neighbors if not atom.name == self.center.name][0]
                c = (x1 - X.coordinates)/norm(x1 - X.coordinates) + X.coordinates
                v1 = X.coordinates - X.neighbors[0].coordinates
                v2 = X.coordinates - X.neighbors[1].coordinates
                a = X.coordinates + (v1 + v2)/norm(v1 + v2)
        return a, b, c
    
    def ep_crds(self, distance, angle):
        angle *= np.pi/180.
        distance *= units.angst_to_au
        ep_global = np.zeros((self.num_ep, 3))
        ep_local = np.zeros((self.num_ep, 3))
        tetas = [0. + self.offset, np.pi + self.offset][:self.num_ep]
        for i, teta in enumerate(tetas):
            # in local coordinate system:
            ep_local[i][0] = np.sin(angle)*np.sin(teta)*distance
            ep_local[i][1] = np.sin(angle)*np.cos(teta)*distance 
            ep_local[i][2] = np.cos(angle)*distance
            # in global
            ep_global[i] = np.dot(self.local_system, ep_local[i])
            ep_global[i] += self.center.coordinates
        logger.debug("Coordinates of extra points:\n%r" % ep_global)
        return ep_global

    def ep_crds_fast(self, distance, angle):
        return fast.ep_crds(distance, angle, self.tetas, self.local_system, self.center.coordinates, self.num_ep)

