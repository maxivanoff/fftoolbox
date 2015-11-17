
# coding: utf-8

# In[1]:

import matplotlib.pyplot as plt
import seaborn as sns
import logging, sys

logger = logging.getLogger(__name__)
lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
logging.basicConfig(filename='./log', filemode='w', level=logging.DEBUG, format=lf)


# In[2]:

import fftoolbox as fftb


# In[3]:

molecule_name = 'methanol'
theory = 'mp2_augccpvtz'
max_rank = 1
radius = 1.0
ranks = range(max_rank+1)


# In[4]:

def get_reports(ChargeModel, grid):
    reports = fftb.Reports()
    for rank in ranks:
        data = {
                'name': molecule_name,
                'theory': theory,
                'sphere params': (rank, radius),
                }
        # create Lebedev Charge Model for the molecule 
        parser = fftb.GDMA(data=data)
        data.update(parser.data.copy())
        molecule = ChargeModel(data)
        # Link charges with the grid
        charges = fftb.LeastSquaresCharges(molecule, grid)
        charges.sites_to_solution()
        reports.add(charges)
    return reports


# ### Set up van der Waals grids

# In[5]:

data = {
        'name': molecule_name,
        'theory': theory,
        'density': 1.5,
        'exclude': ['<xy']
        }
parser = fftb.GaussianCube(data=data)
data.update(parser.data.copy())
grids = fftb.vdwGrids(data)
for grid_name, grid in grids.items():
    grid.build_LEM(filename='%s.pymol' % grid_name)


# ### Atom-centered charges

# In[6]:

molecule = fftb.HM(data)
charges = fftb.LeastSquaresCharges(grid=grids['full vdw'], molecule=molecule)
charges.solve()
reportAC = fftb.Report(name='atom-centered %s' % molecule_name, charges=charges)
print reportAC


# In[7]:

reports = {}
for grid_name, grid in grids.items():
    reports_per_grid = {
        'LM': get_reports(fftb.LM, grid),
        'DLM': get_reports(fftb.DLM, grid),
        }
    reports[grid_name] = reports_per_grid


# In[8]:

for grid_name, grid in grids.items():
    reports_per_grid = reports[grid_name]
    for i, key in enumerate(reports_per_grid['LM'].keys):
        fig = plt.figure(i)
        
        plt.subplot(1,2,1)
        axes = plt.gca()
        plt.title('LM %s' % key)
        #plt.plot(ranks, [ACreport.data[key]]*len(ranks), '-',color='black')
        plt.plot(ranks, reports_per_grid['LM'].data[key],label=grid_name)
        plt.plot(ranks, reports_per_grid['LM'].data[key],'o', color='grey', ms=6)
        ylim = axes.get_ylim()
        plt.legend(loc='upper right')        
        
        plt.subplot(1,2,2)
        axes = plt.gca()
        axes.set_ylim(ylim)
        plt.title('DLM %s' % key)
        plt.plot(ranks, reports_per_grid['DLM'].data[key],label=grid_name)
        plt.plot(ranks, reports_per_grid['DLM'].data[key],'o', color='grey', ms=6)
        plt.legend(loc='upper right')


# In[ ]:



