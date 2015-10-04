import matplotlib.pyplot as plt
import matplotlib
import numpy as np

def table(t, name, sites_names, singulars):
    font = {'size'   : 32}
    matplotlib.rc('font', **font)
    t = np.around(t, decimals=2)
    dc = 0.2 # color depth
    left, width = 0.1, 0.9
    bottom, height = 0.2, 0.9


    M_table = [left, bottom , width, height]
    axM = plt.axes(M_table, frameon =False, xticks=[], yticks=[])
    normal = plt.Normalize(t.min()-dc, t.max()+dc)
    axM.table(cellText=t, loc='center', colWidths = [0.05]*t.shape[1],\
            rowLabels=[s[0] for s in sites_names], colLabels=['%.2f' % s for s in singulars],\
                cellColours=plt.cm.coolwarm(normal(t)))
    plt.savefig('%s.pdf' % (name))
    plt.close()

