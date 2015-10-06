from molecule import Molecule
from parser import GaussianCube, QChem, Gaussian, Mol2
import logging, sys

logger = logging.getLogger(__name__)
lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
logging.basicConfig(level=logging.DEBUG, format=lf)

from xml.etree.ElementTree import Element, SubElement, Comment
from xml.etree import ElementTree
from xml.dom import minidom

def prettify(elem):
    """Return a pretty-printed XML string for the Element.
    """
    rough_string = ElementTree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")

LJ = {
        'CT': (1.9080, 0.1094), # methyl carbon
        'C' : (1.9080, 0.0860), # sp2
        'HC': (1.4870, 0.0157), # connected to CT / methyl hydrogen,
        'H' : (0.6000, 0.0157), # connected to N
        'HO': (0.0000, 0.0000), # connected to OH / hydroxyl
        'N' : (1.8240, 0.1700), # sp2
        'N3': (1.8750, 0.1700), # sp3
        'OH': (1.7210, 0.2104),
        'O' : (1.6612, 0.2100), # C=O
        'O2': (1.6612, 0.2100), # CO2
        'S' : (2.0000, 0.2500),
        }

amber = {
        'C0_CH3': 'CT',
        'H0_CH3_0': 'HC',
        'H0_CH3_1': 'HC',
        'H0_CH3_2': 'HC',
        'N': 'N',
        'S': 'S',
        'O': 'O',
        }

data = {
        'name': 'trans-mesno',
        'theory': 'pbe_def2svp',
        }

# load molecule

fname = 'truncate-2/%s_%s.mol2' % (data['name'], data['theory'])
parser = Mol2(filename=fname)
data.update(parser.data.copy())
molecule = Molecule(data=data)

# write to xml

top = Element('forcefield', name=data['name'])
comment = Comment('Generated from %s' % fname)
top.append(comment)

for site in molecule.sites_noneq:
    siteElem = SubElement(top, 'site', name=site.name)
    try:
        r0, e = LJ[amber[site.name]]
    except KeyError:
        r0, e = 0., 0.
    SubElement(siteElem, 'charge').text = str(site.charge)
    SubElement(siteElem, 'r0').text = str(r0)
    SubElement(siteElem, 'epsilon').text = str(e)

s = prettify(top)
xmlfilename = 'data/forcefields/%s.xml' % fname.split('.')[0]
f = open(xmlfilename, 'w') 
f.write(s)
f.close()

