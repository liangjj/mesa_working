*deck pchdoc
      subroutine pchdoc
c***begin prologue  pchdoc
c***purpose  documentation for pchip, a fortran package for piecewise
c            cubic hermite interpolation of data.
c***library   slatec (pchip)
c***category  e1a, z
c***type      all (pchdoc-a)
c***keywords  cubic hermite interpolation, documentation,
c             monotone interpolation, pchip,
c             piecewise cubic interpolation
c***author  fritsch, f. n., (llnl)
c             lawrence livermore national laboratory
c             p.o. box 808  (l-316)
c             livermore, ca  94550
c             fts 532-4275, (510) 422-4275
c***description
c
c            pchip:  piecewise cubic hermite interpolation package
c
c      this document describes the contents of pchip, which is a
c   fortran package for piecewise cubic hermite interpolation of data.
c   it features software to produce a monotone and "visually pleasing"
c   interpolant to monotone data.  as is demonstrated in reference 4,
c   such an interpolant may be more reasonable than a cubic spline if
c   the data contains both "steep" and "flat" sections.  interpola-
c   tion of cumulative probability distribution functions is another
c   application.  (see references 2-4 for examples.)
c
c
c      all piecewise cubic functions in pchip are represented in
c   cubic hermite form; that is, f(x) is determined by its values
c   f(i) and derivatives d(i) at the breakpoints x(i), i=1(1)n.
c   throughout the package a pch function is represented by the
c   five variables  n, x, f, d, incfd:
c     n     - number of data points;
c     x     - abscissa values for the data points;
c     f     - ordinates (function values) for the data points;
c     d     - slopes (derivative values) at the data points;
c     incfd - increment between successive elements in the f- and
c             d-arrays (more on this later).
c   these appear together and in the same order in all calls.
c
c      the double precision equivalents of the pchip routines are
c   obtained from the single precision names by prefixing the
c   single precision names with a d.  for example, the double
c   precision equivalent of pchim is dpchim.
c
c      the contents of the package are as follows:
c
c   1. determine derivative values.
c
c      note:  these routines provide alternate ways of determining d
c             if these values are not already known.
c
c         pchim -- piecewise cubic hermite interpolation to monotone
c               data.
c               used if the data are monotonic or if the user wants
c               to guarantee that the interpolant stays within the
c               limits of the data.  (see reference 3.)
c
c         pchic -- piecewise cubic hermite interpolation coefficients.
c               used if neither of the above conditions holds, or if
c               the user wishes control over boundary derivatives.
c               will generally reproduce monotonicity on subintervals
c               over which the data are monotonic.
c
c         pchsp -- piecewise cubic hermite spline.
c               produces a cubic spline interpolator in cubic hermite
c               form.  provided primarily for easy comparison of the
c               spline with other piecewise cubic interpolants.  (a
c               modified version of de boor's cubspl, reference 1.)
c
c   2. evaluate, differentiate, or integrate resulting pch function.
c
c      note:  if derivative values are available from some other
c             source, these routines can be used without calling
c             any of the previous routines.
c
c         chfev -- cubic hermite function evaluator.
c               evaluates a single cubic hermite function at an array
c               of points.  used when the interval is known, as in
c               graphing applications.  called by pchfe.
c
c         pchfe -- piecewise cubic hermite function evaluator.
c               used when the interval is unknown or the evaluation
c               array spans more than one data interval.
c
c         chfdv -- cubic hermite function and derivative evaluator.
c               evaluates a single cubic hermite function and its
c               first derivative at an array of points.  used when
c               the interval is known, as in graphing applications.
c               called by pchfd.
c
c         pchfd -- piecewise cubic hermite function and derivative
c               evaluator.
c               used when the interval is unknown or the evaluation
c               array spans more than one data interval.
c
c         pchid -- piecewise cubic hermite integrator, data limits.
c               computes the definite integral of a piecewise cubic
c               hermite function when the integration limits are data
c               points.
c
c         pchia -- piecewise cubic hermite integrator, arbitrary limits.
c               computes the definite integral of a piecewise cubic
c               hermite function over an arbitrary finite interval.
c
c   3. utility routines.
c
c         pchbs -- piecewise cubic hermite to b-spline converter.
c               converts a pch function to b-representation, so that
c               it can be used with other elements of the b-spline
c               package (see bspdoc).
c
c         pchcm -- piecewise cubic hermite, check monotonicity of.
c               checks the monotonicity of an arbitrary pch function.
c               might be used with pchsp to build a polyalgorithm for
c               piecewise c-2 interpolation.
c
c   4. internal routines.
c
c         chfie -- cubic hermite function integral evaluator.
c               (real function called by pchia.)
c
c         chfcm -- cubic hermite function, check monotonicity of.
c               (integer function called by pchcm.)
c
c         pchce -- pchic end derivative setter.
c               (called by pchic.)
c
c         pchci -- pchic initial derivative setter.
c               (called by pchic.)
c
c         pchcs -- pchic monotonicity switch derivative setter.
c               (called by pchic.)
c
c         pchdf -- pchip finite difference formula.
c               (real function called by pchce and pchsp.)
c
c         pchst -- pchip sign testing routine.
c               (real function called by various pchip routines.)
c
c         pchsw -- pchcs switch excursion adjuster.
c               (called by pchcs.)
c
c   the calling sequences for these routines are described in the
c   prologues of the respective routines.
c
c
c      incfd, the increment between successive elements in the f-
c   and d-arrays is included in the representation of a pch function
c   in this package to facilitate two-dimensional applications.  for
c   "normal" usage incfd=1, and f and d are one-dimensional arrays.
c   one would call pchxx (where "xx" is "im", "ic", or "sp") with
c
c              n, x, f, d, 1  .
c
c   suppose, however, that one has data on a rectangular mesh,
c
c         f2d(i,j) = value at (x(i), y(j)),  i=1(1)nx,
c                                            j=1(1)ny.
c   assume the following dimensions:
c
c         real  x(nxmax), y(nymax)
c         real  f2d(nxmax,nymax), fx(nxmax,nymax), fy(nxmax,nymax)
c
c   where  2.le.nx.le.nxmax and 2.le.ny.le.nymax .  to interpolate
c   in x along the line  y = y(j), call pchxx with
c
c              nx, x, f2d(1,j), fx(1,j), 1  .
c
c   to interpolate along the line x = x(i), call pchxx with
c
c              ny, y, f2d(i,1), fy(i,1), nxmax  .
c
c   (this example assumes the usual columnwise storage of 2-d arrays
c    in fortran.)
c
c***references  1. carl de boor, a practical guide to splines, springer-
c                 verlag, new york, 1978 (esp. chapter iv, pp.49-62).
c               2. f. n. fritsch, piecewise cubic hermite interpolation
c                 package, report ucrl-87285, lawrence livermore natio-
c                 nal laboratory, july 1982.  [poster presented at the
c                 siam 30th anniversary meeting, 19-23 july 1982.]
c               3. f. n. fritsch and j. butland, a method for construc-
c                 ting local monotone piecewise cubic interpolants, siam
c                 journal on scientific and statistical computing 5, 2
c                 (june 1984), pp. 300-304.
c               4. f. n. fritsch and r. e. carlson, monotone piecewise
c                 cubic interpolation, siam journal on numerical ana-
c                 lysis 17, 2 (april 1980), pp. 238-246.
c***routines called  (none)
c***revision history  (yymmdd)
c   811106  date written
c   870930  updated reference 3.
c   890414  changed pchmc and chfmc to pchcm and chfcm, respectively,
c           and augmented description of pchcm.
c   891214  prologue converted to version 4.0 format.  (bab)
c   910826  1. revised purpose, clarified role of argument incfd,
c              corrected error in example, and removed redundant
c              reference list.
c           2. added description of pchbs.  (fnf)
c   920429  revised format and order of references.  (wrb,fnf)
c   930505  changed chfiv to chfie.  (fnf)
c***end prologue  pchdoc
c-----------------------------------------------------------------------
c     this is a dummy subroutine, and should never be called.
c
c***first executable statement  pchdoc
      return
c------------- last line of pchdoc follows -----------------------------
      end
