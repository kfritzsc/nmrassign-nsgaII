# !/usr/bin/env python
# setup file for nmrassign

from distutils.core import setup
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# get long description from README
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='nmrassign',
      version='0.1-dev',
      description='Python tools for assigning protein NMR data',
      long_description=long_description,
      author='Keith Fritzsching',
      author_email='kfritzsc@brandeis.edu',
      url='https://github.com/kfritzsc/nmrassign',
      packages=['nmrassign'],
      requires=['numpy'],
      license='New BSD License',
      classifiers=['Intended Audience :: Science/Research',
                   'Intended Audience :: Developers',
                   'License :: OSI Approved :: BSD License',
                   'Programming Language :: Python :: 2',
                   'Programming Language :: Python :: 2.7',
                   'Programming Language :: Python :: 3',
                   'Programming Language :: Python :: 3.5',
                   'Topic :: Scientific/Engineering',
                   'Operating System :: MacOS :: MacOS X',
                   'Operating System :: Microsoft :: Windows',
                   'Operating System :: POSIX :: Linux'],)
