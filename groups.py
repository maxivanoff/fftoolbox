from multipole import Multipole
import numpy as np
import logging
logger = logging.getLogger(__name__)

class CarbonylGroup(Multipole):

    def __init__(self, center, count, name=None, ff=None, sym=False):
        self.carbon = center
        self.oxygen = filter(lambda a: a.element=='O', center.neighbors)[0]
        self.hydrogens = filter(lambda a: a.element=='H', center.neighbors)
        Multipole.__init__(self, name, ff)
        self.atoms = [self.carbon, self.oxygen]
        self.name = 'CO_%i' % count
        #self.carbon.name = 'C%i_CO' % count 
        #self.oxygen.name = 'O%i_CO' % count 
        self.oxygen.name = 'O'
        if len(self.hydrogens)==1:
            self.hydrogens[0].name='H%i_CO' % count


class BuriedGroup(Multipole):

    def __init__(self, center, count, name=None, ff=None, sym=False):
        self.center = center
        name = 'buried %s' % center.element
        self.nonhydrogens = filter(lambda a: a.element!='H', center.neighbors)
        self.hydrogens = filter(lambda a: a.element=='H', center.neighbors)
        Multipole.__init__(self, name, origin=self.get_origin())
        self.atoms = [center] + self.hydrogens
        cname = self.center.element
        if self.hydrogens:
            neighbor_name = 'H'
            neighbors = self.hydrogens
        else:
            neighbor_name = self.nonhydrogens[0].element
            neighbors = self.nonhydrogens
        if [neighbor_name]*3==[a.element for a in neighbors]:
            self.name = cname + neighbor_name + '3_%i' % (count)
            self.type = cname + neighbor_name + '3'
            self.center.name = '%s%i_%s%s3' % (cname, count, cname, neighbor_name)
        elif [neighbor_name]*2==[a.element for a in neighbors]:
            self.name = cname + neighbor_name + '2_%i' % (count)
            self.type = cname + neighbor_name + '2'
            self.center.name = '%s%i_%s%s2'% (cname, count, cname, neighbor_name)
        elif [neighbor_name]*1==[a.element for a in neighbors]:
            self.name = cname + neighbor_name + '1_%i' % (count)
            self.type = cname + neighbor_name + '1'
            self.center.name = '%s%i_%s%s1'% (cname, count, cname, neighbor_name)
        elif [neighbor_name]*4 == [a.element for a in neighbors]:
            self.name = cname + neighbor_name+'4'
            self.type = cname + neighbor_name+'4'
            self.center.name = cname
        else:
            raise ValueError('Not a buired group')
        for i,a in enumerate(neighbors):
            if sym is True:
                a.name = neighbor_name + center.name[1:]
            else:
                a.name = neighbor_name + center.name[1:] + '_' + str(i)
        logger.debug('Buried atom created with center @ %s and %i hydrogen atoms: %r' %(self.center.name, len(self.hydrogens), [a.name for a in self.hydrogens]))

    def get_origin(self, alpha=None, pos='non-hydrogen'):
        if alpha:
            # origin at virtual site
            return self.get_CX(alpha) + self.center.coordinates
        if pos=='non-hydrogen':
            # origin at non hydrogen atom
            atoms = self.nonhydrogens
        if pos=='hydrogen':
            atoms = self.hydrogens

        origin = np.zeros(3)
        for a in atoms:
            origin += a.coordinates
        if not len(atoms) == 0:
            return origin/len(atoms)
        else:
            return origin

    def get_CX(self, alpha):
        CH = np.zeros(3)
        for h in self.hydrogens:
            CH += self.center.coordinates - h.coordinates
        return -CH/(len(self.hydrogens)-abs(alpha))

    def get_dipole_norm(self):
        dipole = Multipole.get_dipole(self)
        CX = np.zeros(3)
        for a in self.nonhydrogens:
            CX += a.coordinates - self.center.coordinates
        sign = np.sign(np.dot(dipole, CX))
        return sign*np.linalg.norm(dipole)

