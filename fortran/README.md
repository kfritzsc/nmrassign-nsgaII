Combination NSGA-II/MC NMR Assignment
=====================================

Algorithm for automatic SSNMR resonance assignment

With gcc-gfortrain installed, type make to compile the Fortran code.

How to run
----------

When executing this program, first the terminal will display: 'Enter control
file name:'. This is where you input the name and path of the control file.
Then you can verify the input information included in the control file at the
command window.  Wait for sometime
(normally 5 - 60 minutes) and check the output files (they will be at the
file path you set).

Input files
-----------

Prepare the input files before using this program (Examples of the input files
are in the "tests" folder):
1. Amino acid sequence
2. Peak lists
3. Connection table
4. Control file (the program will upload the input files and adjust the
   parameters according to this file) The control file is explained in detail
   below.

Output files
------------

The NSGA-II/MC program will output the results in two different ways:
1. output file(data)	- the M results (M is the group size) with the peak
number are shown in this file. This file is not for the users to read, but can
be used as the initial data file (see instructions for the control file) and
for the post-processing MATLAB program.
2. output files(table)	- each result is outputted to an independent file with
the chemical shifts of each residue shown.

Control file
------------

- sequence file: the name (and path) of the sequence file
- number of spectra: if the number is N, then there should be N columns
  showing the names of the spectra below.
- spectra 1 & 2 & more: the names (and paths) of the spectra
- connection table: the name (and path) of the connection table file
- initial data file: for initializing the individuals. If there are no initial
  data ("NULL" as input), then the program will generate the initial data
  randomly.
- output file(data): the name (and path) of the output data
- output files(table): the name (and path) of the output tables
- group size: the number of the results (larger numbers will increase the
  computation time, normally smaller than 100)
- gene pool size: the number of individuals that are chosen in the "gene
  pool", should be smaller than the group size
- number of steps:  the number of steps for "lucky ratio" variation
- number of NSGA-II attempts:  the number of NSGA-II loops in each step
- number of MC attempts: the number of Monte-Carlo loops in each step
- number of free steps: the number of steps that the algorithm is running
  using "unrestricted Pareto-order strategy"
- mutation rate : the rate of the mutation operator for the genetic algorithm
- additional mutation rate: the rate of an additional mutation operator for
  the genetic algorithm
- crossover rate: the rate of the crossover operator for the genetic algorithm
- null probability: the rate to change one assignment to zero (i.e. null
  assignment)
- minimum w1-w4:  minimum weighting factors w1-w4 for Monte-Carlo attempts
- maximum w1-w4:  maximum weighting factors w1-w4 for Monte-Carlo attempts
