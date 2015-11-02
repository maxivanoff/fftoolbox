import numpy as np
import logging

from multipole import Multipole

__author__ = "Maxim Ivanov"
__email__ = "maxim.ivanov@marquette.edu"

logger = logging.getLogger(__name__)

class CarbonOxygenGroup(object):

    def __init__(self, center, count, name=None, sym=False):
        if count == 0: count = ''
        oxygens = filter(lambda a: a.element=='O', center.neighbors)
        if len(oxygens) == 1:
            self.name = 'CO'
        if len(oxygens) == 2:
            self.name = 'CO2'
        name = 'C%s_%s' % (count, self.name)
        center.set_name(name)
        for i, oxygen in enumerate(oxygens):
            if sym is False:
                name = 'O%s_%s-%i' % (count, self.name, i)
                oxygen.set_name(name)
            else:
                name = 'O%s_%s' % (count, self.name)
                oxygen.set_name(name)
        #self.atoms = oxygens + [center]
        #Multipole.__init__(self, name=self.name)
        #self.hydrogens = filter(lambda a: a.element=='H', center.neighbors)
        logger.debug('%s group identified with center %s, %i oxygens' % (self.name, center.name, len(oxygens)))

class BuriedGroup(object):
    """
    Methyl or methylene group
    """

    def __init__(self, center, count, name=None, sym=False):
        if count == 0: count = ''
        if not center.element == 'C':
            raise ValueError('Should be a buried carbon group')
        hydrogens = filter(lambda a: a.element=='H', center.neighbors)
        num_h = len(hydrogens)
        self.name = 'CH%i' % num_h
        name = 'C%s_CH%i' % (count, num_h)
        center.set_name(name)
        # set hydrogens names
        for i, a in enumerate(hydrogens):
            if sym is True:
                name = 'H%s' % center.name[1:]
                a.set_name(name)
            else:
                name = 'H%s-%i' % (center.name[1:], i)
                a.set_name(name)
        self.atoms = [center] + hydrogens
        logger.debug('Buried atom created with center @ %s and %i hydrogen atoms: %r' %(center.name, len(hydrogens), [a.name for a in hydrogens]))

