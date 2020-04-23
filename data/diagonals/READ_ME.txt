This folder contains, for each of the 90 instances of Mean-Variance (MV) portfolio optimization problem with minimum and maximum buy-in thresholds available at

  http://www.di.unipi.it/optimize/Data/MV.html

some diagonal matrices D with the property that Q - D is positive semidefinite, where Q is the variance-covariance matrix of the assets (i.e., the Hessian of the MV problem). This is useful to apply Perspective Relaxation techniques: the quality of the bound is tied to the "quality" of the matrix.

Two different ways to compute the matrices are proposed in the articles

S) A. Frangioni, C. Gentile "SDP Diagonalizations and Perspective Cuts for a Class of Nonseparable MIQP" Operations Research Letters, 35(2):181–185, 2007

L) X. Zheng, X. Sun, D. Li. "Improving the Performance of MIQP Solvers for Quadratic Programs with Cardinality and Minimum Threshold Constraints: A Semidefinite Program Approach" INFORMS Journal on Computing, 26(4):690–703, 2014

Both require solving a SemiDefinite Program (SDP): that in [S] is "smaller", thus the matrices are indicated with D_s in [L], where a "larger" SDP is proposed which produces matrices D_l that yield a better root node bound. Sometimes, using a convex combination of the two matrices results in better overall B&C performances; in [L], the matrices

   D_c = 0.5 * D_l + 0.5 * D_s

are also tested. It is of course useless to separately distribute them.

While the SDP in [S] only uses the Hessian Q, that in [L] uses all the constraints of the problem. Hence, matrices D_l depend, for the same instance, also on the bound on the maximum size of the portfolio, k. The diagonals are stored in files, with the same name as the oringinal MV instance file name, in the subdirectories "s", "l-10" and "l-n". "s" contains the D_s diagonals (independent from k), "l-10" contains the D_l for k = 10, and "l-n" contains the D_l for k = number of assets, i.e., for (MV) without cardinality constraints. Note that the (MV) instances also have "min buy-in constraints" imposing a minimum amount of each asset to be bought; these imply that a maximum of 20 assets can be bought, even if no explicit cardinality constraint is imposed.

The format of the files is very simple:

< n >
for i = 1 to n
 < D[ i ][ i ] = diagonal entry of asset i >

each number on a separate line.

These matrices are used for computational tests in

A. Frangioni, F. Furini, C. Gentile "Improving the Approximated Projected Perspective Reformulation by Dual Information" Technical Report, Dipartimento di Informatica, Universita di Pisa, 2016



