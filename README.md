nmrassign
=========
Fortran and Python tools for assigning protein NMR data with multi-objective optimization.

Tools to help assign protein (solid-state) NMR spectra with resonance
overlap and missing peaks.  The project is based on the methods of Tycko *et.
al.* [1-2] and Yang *et. al.* [3].

The primary motivation for the project is to ease the production of input data
for the multi-objective optimization routines introduced by Yang *et. al.*.

Install
-------

1. **Install python tools.**
   ```
   python setup.py install
   ```

2. **Compile Fortran tools.**
   ```
   cd <project path>/fortran/nsga2/
   make
   ```
   The Makefile works well with gcc-gfortran from brew.

References
----------

1. Tycko, R. and Hu, K. N. *J. Magn. Reson.* 205, 304-314 (2010).
2. Hu, K. N.; Qiang, W.; Tycko, R. *J. Biomol. NMR.* 50, 267-276 (2011).
3. Yang, Y.; Fritzsching, K. J.; Hong, M. *J. Biomol. NMR.* 57, 281-296
   (2013).
