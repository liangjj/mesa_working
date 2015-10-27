 input parameter bounds flag vector nbdk(k) =
 2 2 2
 bnds 36, nbdsum =  4, par bd vec, after possible modification by par flag vec, =
 0 2 2
 iters  =  40 number of allowed iterations
 relphi =  1.00E-13, convergence when relative drop in phi less than input relphi
 phimin =  1.00E-25, convergence when phi less than nonzero phimin
 grdzmin=  1.00E-20, convergence when grad.z between 0.0' and -grdzmin
zsqmin =  1.00E-24exit when squared mag of z vector less than zsqmin
1problem  0 begins here
 parameter flag vec nuflag(k), 0 (1) to vary (hold fixed)' parameter v(k)
 parameter bound flag vec nbd(k), 0 (v(k) unbounded), 1 (v(k) has lower bound only), 2 (v(k) bounded below and above)
 3' (v(k) has upper bound only)
 parameter lower bound vec bl(k), format(6e12.6),  blank component when v(k) is unbounded below or held fixed
 parameter upper bound vec bu(k), format(6e12.6),  blank component when v(k) is unbounded above or held fixed
 opt 180, nuflag(k) nbd(k)    bl(k)    parameter vector v(k)       bu(k)        k
                1        0   1.000000E+00   1.00000000000000E+00   1.000000E+00   1
                0        2   3.415000E+00   3.42500000000000E+00   3.430000E+00   2
                0        2   2.000000E-02   2.55000000000000E-02   3.000000E-02   3
 ****    ****    ****    ****    ****    ****    ******
 small z vector exit, squared mag z vec =  6.83E-29 is less than input zsqmin =  1.00E-24
 ****    ****    ****    ****    ****    ****    ******
 230, phi =    3.40215395594669E-12, parameter vector u =
      , end of iteration 11

 statistical output, standard deviations
    2.31E-08  -5.98E-09
    6.38E-08   3.95E-08
 opt 230, parameter vec v =
    1.00000000000000E+00    3.42639028058237E+00    2.55490020164543E-02
 230, gradient squared mag =  6.87E-13,gradient vec =
  8.21E-07 -1.11E-07
 problem 0, running time =  0.000000E+00
