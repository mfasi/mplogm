# Multiprecision matrix logarithm

This repository contains the code that produces the figures and tables in the
technical report:

M. Fasi and N. J. Higham, Multiprecision Algorithms for Computing the Matrix
Logarithm. Technical Report 2017.16, Manchester Institute for Mathematical
Sciences, The University of Manchester, UK, May 2017.

## Dependencies

The code in this repository requires the Advanpix Multiprecision Computing
Toolbox for MATLAB (www.advanpix.com). For inclusion into LaTeX docuements, the
figures produced by the MATLAB can be exported in TikZ format using the script
`matlab2tikz` (http://uk.mathworks.com/matlabcentral/fileexchange/22022) and
setting to `true` the variable `tikz_out` in `run_test.m`.

## Execution of the experiments

Assuming that the system is configured correctly, the command `run_tests` will
run the all suite of numerical experiments.