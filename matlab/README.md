postproc_assign.m
=================
Post-processing code for assignment data from NSGA2-Assign and Combo 
NSGA2/ MC-Assign Programs.

This version is created using MATLAB R2009b

Description
-----------

This program collects the Pareto-order-1 solutions from the outputs of
NSGA2-Assign and Combo NSGA2/MC-Assign programs, and calculate the 
assignment probabilities according to these solutions. The assignment 
results with larger than 90% probability will be organized and output 
to a text file.
