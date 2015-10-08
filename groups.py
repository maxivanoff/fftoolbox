import numpy as np
import logging

from multipole import Multipole

__author__ = "Maxim Ivanov"
__email__ = "maxim.ivanov@marquette.edu"

logger = logging.getLogger(__name__)

class CarbonOxygenGroup(Multipole):

    def __init__(self, center, count, name=None, sym=False):
        if count == 0: count = ''
        self.carbon = center
        self.oxygens = filter(lambda a: a.element=='O', center.neighbors)
        if len(self.oxygens) == 1:
            self.name = 'CO'
        if len(self.oxygens) == 2:
            self.name = 'CO2'
        self.carbon.name = 'C%s_%s' % (count, self.name)
        for i, oxygen in enumerate(self.oxygens):
            if sym is False:
                oxygen.name = 'O%s_%s-%i' % (count, self.name, i)
            else:
                oxygen.name = 'O%s_%s' % (count, self.name)
        self.atoms = self.oxygens + [self.carbon]
        Multipole.__init__(self, name=self.name)
        self.hydrogens = filter(lambda a: a.element=='H', center.neighbors)
        logger.debug('%s group identified with carbon %s, %i oxygens and %i hydrogens' % (self.name, self.carbon.name, len(self.oxygens), len(self.hydrogens)))


class BuriedGroup(Multipole):

    def __init__(self, center, count, name=None, sym=False):
        if count == 0: count = ''
        self.center = center
        self.nonhydrogens = filter(lambda a: not a.element=='H', center.neighbors)
        self.hydrogens = filter(lambda a: a.element=='H', center.neighbors)
        self.atoms = [center] + self.hydrogens
        cname = self.center.element
        if self.hydrogens:
            neighbor_name = 'H'
            neighbors = self.hydrogens
        else:
            neighbor_name = self.nonhydrogens[0].element
            neighbors = self.nonhydrogens
        if [neighbor_name]*3==[a.element for a in neighbors]:
            self.name = cname + neighbor_name + '3'
            self.center.name = '%s%s_%s%s3' % (cname, count, cname, neighbor_name)
        elif [neighbor_name]*2==[a.element for a in neighbors]:
            self.name = cname + neighbor_name + '2'
            self.center.name = '%s%s_%s%s2'% (cname, count, cname, neighbor_name)
        elif [neighbor_name]*1==[a.element for a in neighbors]:
            self.name = cname + neighbor_name + '1'
            self.center.name = '%s%s_%s%s1'% (cname, count, cname, neighbor_name)
        elif [neighbor_name]*4 == [a.element for a in neighbors]:
            self.name = cname + neighbor_name+'4'
            self.center.name = '%s%s_%s%s4'% (cname, count, cname, neighbor_name)
        else:
            raise ValueError('Not a buired group')
        for i,a in enumerate(neighbors):
            if sym is True:
                a.name = neighbor_name + center.name[1:]
            else:
                a.name = neighbor_name + center.name[1:] + '-' + str(i)
        Multipole.__init__(self, name=self.name)
        logger.debug('Buried atom created with center @ %s and %i hydrogen atoms: %r' %(self.center.name, len(self.hydrogens), [a.name for a in self.hydrogens]))

