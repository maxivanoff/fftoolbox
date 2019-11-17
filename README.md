#### *fftoolbox*: Python library to optimize force fields

* fits atom-centered point charges to the reference electrostatic potential
* reference grid has a form of thick van der Waals surface (Gaussian .cub file required) 
* computes multipole moments for the derived charges (in spherical or cartesian representation)
* written in Python 2.7, conversion to Python 3.0 is in progress
* requires Gaussian log files in data to run

Matrix of inversed distances is computed using Cython (see fast.pyx)

run `python setup.py build_ext --inplace` to compile fast.pyx

run `python run-simple.py` for a simple example

works that use *fftoolbox*:

*J. Phys. Chem. A*, 119 (8), 1422 (2015)

*J. Chem. Phys.*, 143, 134102 (2015)

