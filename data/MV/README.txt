This is a distribution of instances of the Mean-Variance problem with minimum buy-in constraints,
i.e., the following Mixed-Integer Quadratic program

min x^T Q x

    \sum_i^n \mu_i x_i >= \rho

    \sum_i^n x_i = 1

    l_i y_i <= x_i <= u_i y_i		i = 1, ..., n

    y_i \in {0, 1}			i = 1, ..., n

where Q is a Positive Semi-Definite n \times n matrix, intended to be the variance-covariance matrix
of n marketable assets. Each \mu_i is the expected return of the corresponding assets during the investment
period, and the problem seeks to define the portfolio which minimizes risk (as measured by the covariance of
the portfolio) while providing the required level \rho of return. This would be a convex quadratic problem (= easy)
if it were not for the "complicating" constraint that each asset i has a nonzero minimum quantity l_i of that asset
that can be bought, if the asset is bought at all; there also is a maximum quantity u_i, but this (in itself) does
not create particular problems.

--------------------------------------------------------------
The data of each instance "*" is stored in 4 files:

*.txt

contains:
a) number of assets (n)
b) n rows with two values of type double:
   - the first value is \mu_i (the return value of asset i)
   - the second value must be ignored

*.rho

contains the value of \rho

*.bds

contains n rows with two entries:
- the first entry is l_i
- the second entry is u_i

*.mat

contains the value n (this is a redundant information) and then the matrix Q, that must be symmetric and positive semi-definite, in row-wise format; note that *all* the matrix is there, and not only a triangular part (which would be enough since the matrix is symmetric).

-----------------------------------------------------------------------

The present distribution contains the following files and directories:

README.txt

this file.

BestUBLB.txt

Contains the best upper bound (value of a feasible solution) and lower bound known to us on the optimal value of each instance; since all these instances are solvable by proper methods, these are always less than 0.01% far apart.

size200/
size300/
size400/

are folders containing each 30 different instances of the problem with n = 200, 300 and 400, respectively. Each folder contains three types of instances:

- pard???_?
  instances with diagonal dominance ratio  = 0.63

- orl???_005_?
  instances with diagonal dominance ratio = - 0.05

- orl???_05_?
  instances with diagonal dominance ratio = -0.5

mat_params/

is a folder containing 90 *.param files, one for each instance. These are the parameter files for the random generator described in

P.M. Pardalos and G.P. Rodgers "Computing aspects of a branch and bound
algorithm for quadratic zero-one programming" Computing 144, 45-131 (1990).

and each file contains the parameters that have been used to obtain the Q matrices (*.mat file) associated with the corresponding instances.

params/

is a folder containing the random generator genfin.c, written by Claudio Gentile, that generates the other data of the problem, plus one .param file for each instances containing the corresponding parameters of the generator. These are:

- size of the instance (int)
- average lower bound (double)
- average upper bound (double)
- maximum bound variation (double)
- min expected return (double)
- max expected return (double)
- seed (unsigned int)

The generator produces uniformly distributed data with these parameters.

-------------------------------------------------------------------------

Last updated June 10, 2008
