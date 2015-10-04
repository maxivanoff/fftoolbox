from fftoolbox2.atom import Site
from fftoolbox2.multipole import Multipole
import lebedev_write
import numpy as np
import logging
logger = logging.getLogger(__name__)

class LebedevMolecule(Multipole):

    def __init__(self, data):
        Multipole.__init__(self, data['name'], multipoles=data['multipoles'])
        self.num_sites = data['points']
        self.radius = data['radius']
        self.set_sites()
        self.set_sym_sites()
        self.set_multipole_matrix()
        logger.info("Lebedev molecule is created.\nNumber of charged sites: %i\nRadius: %.1f" % (len(self.sites), self.radius))

    def set_sites(self):
        sites = lebedev_write.Lebedev(self.num_sites)
        self.atoms = []
        self.weights = np.zeros(len(sites))
        for i, site in enumerate(sites):
            xyz = np.array(site[0:3])*self.radius
            w = site[3]*np.sqrt(4*np.pi)
            self.weights[i] = w
            s = Site(name='%i' % i, coordinates=xyz)
            self.atoms.append(s)
        self.atoms = (self.atoms)

            

