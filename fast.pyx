from __future__ import division
import numpy as np
cimport numpy as np
cimport cython

__author__ = "Maxim Ivanov"
__email__ = "maxim.ivanov@marquette.edu"

@cython.boundscheck(False)
def set_inversed(np.ndarray[np.float64_t, ndim=2] grid, np.ndarray[np.float64_t, ndim=2] sites,\
        int n_red_sites, np.ndarray[np.int_t, ndim=1] sym_sites_ind):
    assert grid.dtype==np.float64 and sites.dtype==np.float64 and sym_sites_ind.dtype==np.int

    cdef int n_points = grid.shape[0]
    cdef int n_sites = sites.shape[0]
    cdef int p, s

    cdef double x
    cdef double y
    cdef double z
    
    cdef np.ndarray[np.float64_t, ndim=2] A = np.zeros([n_points, n_sites], dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=2] Ared = np.zeros([n_points, n_red_sites], dtype=np.float64)
    
    for p in range(n_points):
        for s in range(n_sites):
            x = grid[<unsigned int>(p), 0] - sites[<unsigned int>(s), 0]
            y = grid[<unsigned int>(p), 1] - sites[<unsigned int>(s), 1]
            z = grid[<unsigned int>(p), 2] - sites[<unsigned int>(s), 2]
            A[p, s] = 1/np.sqrt(x*x+y*y+z*z)

    for p in range(n_points):
        for s in range(n_sites):
            Ared[p, sym_sites_ind[s]] += A[p,s]

    return Ared

