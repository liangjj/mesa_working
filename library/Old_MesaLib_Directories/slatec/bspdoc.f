*deck bspdoc
      subroutine bspdoc
c***begin prologue  bspdoc
c***purpose  documentation for bspline, a package of subprograms for
c            working with piecewise polynomial functions
c            in b-representation.
c***library   slatec
c***category  e, e1a, k, z
c***type      all (bspdoc-a)
c***keywords  b-spline, documentation, splines
c***author  amos, d. e., (snla)
c***description
c
c     abstract
c         bspdoc is a non-executable, b-spline documentary routine.
c         the narrative describes a b-spline and the routines
c         necessary to manipulate b-splines at a fairly high level.
c         the basic package described herein is that of reference
c         5 with names altered to prevent duplication and conflicts
c         with routines from reference 3.  the call lists used here
c         are also different.  work vectors were added to ensure
c         portability and proper execution in an overlay environ-
c         ment.  these work arrays can be used for other purposes
c         except as noted in bspvn.  while most of the original
c         routines in reference 5 were restricted to orders 20
c         or less, this restriction was removed from all routines
c         except the quadrature routine bsqad.  (see the section
c         below on differentiation and integration for details.)
c
c         the subroutines referenced below are single precision
c         routines.  corresponding double precision versions are also
c         part of the package, and these are referenced by prefixing
c         a d in front of the single precision name.  for example,
c         bvalu and dbvalu are the single and double precision
c         versions for evaluating a b-spline or any of its deriva-
c         tives in the b-representation.
c
c                ****description of b-splines****
c
c     a collection of polynomials of fixed degree k-1 defined on a
c     subdivision (x(i),x(i+1)), i=1,...,m-1 of (a,b) with x(1)=a,
c     x(m)=b is called a b-spline of order k.  if the spline has k-2
c     continuous derivatives on (a,b), then the b-spline is simply
c     called a spline of order k.  each of the m-1 polynomial pieces
c     has k coefficients, making a total of k(m-1) parameters.  this
c     b-spline and its derivatives have m-2 jumps at the subdivision
c     points x(i), i=2,...,m-1.  continuity requirements at these
c     subdivision points add constraints and reduce the number of free
c     parameters.  if a b-spline is continuous at each of the m-2 sub-
c     division points, there are k(m-1)-(m-2) free parameters; if in
c     addition the b-spline has continuous first derivatives, there
c     are k(m-1)-2(m-2) free parameters, etc., until we get to a
c     spline where we have k(m-1)-(k-1)(m-2) = m+k-2 free parameters.
c     thus, the principle is that increasing the continuity of
c     derivatives decreases the number of free parameters and
c     conversely.
c
c     the points at which the polynomials are tied together by the
c     continuity conditions are called knots.  if two knots are
c     allowed to come together at some x(i), then we say that we
c     have a knot of multiplicity 2 there, and the knot values are
c     the x(i) value.  if we reverse the procedure of the first
c     paragraph, we find that adding a knot to increase multiplicity
c     increases the number of free parameters and, according to the
c     principle above, we thereby introduce a discontinuity in what
c     was the highest continuous derivative at that knot.  thus, the
c     number of free parameters is n = nu+k-2 where nu is the sum
c     of multiplicities at the x(i) values with x(1) and x(m) of
c     multiplicity 1 (nu = m if all knots are simple, i.e., for a
c     spline, all knots have multiplicity 1.)  each knot can have a
c     multiplicity of at most k.  a b-spline is commonly written in the
c     b-representation
c
c               y(x) = sum( a(i)*b(i,x), i=1 , n)
c
c     to show the explicit dependence of the spline on the free
c     parameters or coefficients a(i)=bcoef(i) and basis functions
c     b(i,x).  these basis functions are themselves special b-splines
c     which are zero except on (at most) k adjoining intervals where
c     each b(i,x) is positive and, in most cases, hat or bell-
c     shaped.  in order for the nonzero part of b(1,x) to be a spline
c     covering (x(1),x(2)), it is necessary to put k-1 knots to the
c     left of a and similarly for b(n,x) to the right of b.  thus, the
c     total number of knots for this representation is nu+2k-2 = n+k.
c     these knots are carried in an array t(*) dimensioned by at least
c     n+k.  from the construction, a=t(k) and b=t(n+1) and the spline is
c     defined on t(k).le.x.le.t(n+1).  the nonzero part of each basis
c     function lies in the  interval (t(i),t(i+k)).  in many problems
c     where extrapolation beyond a or b is not anticipated, it is common
c     practice to set t(1)=t(2)=...=t(k)=a and t(n+1)=t(n+2)=...=
c     t(n+k)=b.  in summary, since t(k) and t(n+1) as well as
c     interior knots can have multiplicity k, the number of free
c     parameters n = sum of multiplicities - k.  the fact that each
c     b(i,x) function is nonzero over at most k intervals means that
c     for a given x value, there are at most k nonzero terms of the
c     sum.  this leads to banded matrices in linear algebra problems,
c     and references 3 and 6 take advantage of this in con-
c     structing higher level routines to achieve speed and avoid
c     ill-conditioning.
c
c                     ****basic routines****
c
c     the basic routines which most casual users will need are those
c     concerned with direct evaluation of splines or b-splines.
c     since the b-representation, denoted by (t,bcoef,n,k), is
c     preferred because of numerical stability, the knots t(*), the
c     b-spline coefficients bcoef(*), the number of coefficients n,
c     and the order k of the polynomial pieces (of degree k-1) are
c     usually given.  while the knot array runs from t(1) to t(n+k),
c     the b-spline is normally defined on the interval t(k).le.x.le.
c     t(n+1).  to evaluate the b-spline or any of its derivatives
c     on this interval, one can use
c
c                  y = bvalu(t,bcoef,n,k,id,x,inbv,work)
c
c     where id is an integer for the id-th derivative, 0.le.id.le.k-1.
c     id=0 gives the zero-th derivative or b-spline value at x.
c     if x.lt.t(k) or x.gt.t(n+1), whether by mistake or the result
c     of round off accumulation in incrementing x, bvalu gives a
c     diagnostic.  inbv is an initialization parameter which is set
c     to 1 on the first call.  distinct splines require distinct
c     inbv parameters.  work is a scratch vector of length at least
c     3*k.
c
c     when more conventional communication is needed for publication,
c     physical interpretation, etc., the b-spline coefficients can
c     be converted to piecewise polynomial (pp) coefficients.  thus,
c     the breakpoints (distinct knots) xi(*), the number of
c     polynomial pieces lxi, and the (right) derivatives c(*,j) at
c     each breakpoint xi(j) are needed to define the taylor
c     expansion to the right of xi(j) on each interval xi(j).le.
c     x.lt.xi(j+1), j=1,lxi where xi(1)=a and xi(lxi+1)=b.
c     these are obtained from the (t,bcoef,n,k) representation by
c
c                call bsppp(t,bcoef,n,k,ldc,c,xi,lxi,work)
c
c     where ldc.ge.k is the leading dimension of the matrix c and
c     work is a scratch vector of length at least k*(n+3).
c     then the pp-representation (c,xi,lxi,k) of y(x), denoted
c     by y(j,x) on each interval xi(j).le.x.lt.xi(j+1), is
c
c     y(j,x) = sum( c(i,j)*((x-xi(j))**(i-1))/factorial(i-1), i=1,k)
c
c     for j=1,...,lxi.  one must view this conversion from the b-
c     to the pp-representation with some skepticism because the
c     conversion may lose significant digits when the b-spline
c     varies in an almost discontinuous fashion.  to evaluate
c     the b-spline or any of its derivatives using the pp-
c     representation, one uses
c
c                y = ppval(ldc,c,xi,lxi,k,id,x,inppv)
c
c     where id and inppv have the same meaning and usage as id and
c     inbv in bvalu.
c
c     to determine to what extent the conversion process loses
c     digits, compute the relative error abs((y1-y2)/y2) over
c     the x interval with y1 from ppval and y2 from bvalu.  a
c     major reason for considering ppval is that evaluation is
c     much faster than that from bvalu.
c
c     recall that when multiple knots are encountered, jump type
c     discontinuities in the b-spline or its derivatives occur
c     at these knots, and we need to know that bvalu and ppval
c     return right limiting values at these knots except at
c     x=b where left limiting values are returned.  these values
c     are used for the taylor expansions about left end points of
c     breakpoint intervals.  that is, the derivatives c(*,j) are
c     right derivatives.  note also that a computed x value which,
c     mathematically, would be a knot value may differ from the knot
c     by a round off error.  when this happens in evaluating a dis-
c     continuous b-spline or some discontinuous derivative, the
c     value at the knot and the value at x can be radically
c     different.  in this case, setting x to a t or xi value makes
c     the computation precise.  for left limiting values at knots
c     other than x=b, see the prologues to bvalu and other
c     routines.
c
c                     ****interpolation****
c
c     bintk is used to generate b-spline parameters (t,bcoef,n,k)
c     which will interpolate the data by calls to bvalu.  a similar
c     interpolation can also be done for cubic splines using bint4
c     or the code in reference 7.  if the pp-representation is given,
c     one can evaluate this representation at an appropriate number of
c     abscissas to create data then use bintk or bint4 to generate
c     the b-representation.
c
c               ****differentiation and integration****
c
c     derivatives of b-splines are obtained from bvalu or ppval.
c     integrals are obtained from bsqad using the b-representation
c     (t,bcoef,n,k) and ppqad using the pp-representation (c,xi,lxi,
c     k).  more complicated integrals involving the product of a
c     of a function f and some derivative of a b-spline can be
c     evaluated with bfqad or pfqad using the b- or pp- represen-
c     tations respectively.  all quadrature routines, except for ppqad,
c     are limited in accuracy to 18 digits or working precision,
c     whichever is smaller.  ppqad is limited to working precision
c     only.  in addition, the order k for bsqad is limited to 20 or
c     less.  if orders greater than 20 are required, use bfqad with
c     f(x) = 1.
c
c                      ****extrapolation****
c
c     extrapolation outside the interval (a,b) can be accomplished
c     easily by the pp-representation using ppval.  however,
c     caution should be exercised, especially when several knots
c     are located at a or b or when the extrapolation is carried
c     significantly beyond a or b.  on the other hand, direct
c     evaluation with bvalu outside a=t(k).le.x.le.t(n+1)=b
c     produces an error message, and some manipulation of the knots
c     and coefficients are needed to extrapolate with bvalu.  this
c     process is described in reference 6.
c
c                ****curve fitting and smoothing****
c
c     unless one has many accurate data points, direct inter-
c     polation is not recommended for summarizing data.  the
c     results are often not in accordance with intuition since the
c     fitted curve tends to oscillate through the set of points.
c     monotone splines (reference 7) can help curb this undulating
c     tendency but constrained least squares is more likely to give an
c     acceptable fit with fewer parameters.  subroutine fc, des-
c     cribed in reference 6, is recommended for this purpose.  the
c     output from this fitting process is the b-representation.
c
c              **** routines in the b-spline package ****
c
c                      single precision routines
c
c         the subroutines referenced below are single precision
c         routines. corresponding double precision versions are also
c         part of the package and these are referenced by prefixing
c         a d in front of the single precision name. for example,
c         bvalu and dbvalu are the single and double precision
c         versions for evaluating a b-spline or any of its deriva-
c         tives in the b-representation.
c
c     bint4 - interpolates with splines of order 4
c     bintk - interpolates with splines of order k
c     bsqad - integrates the b-representation on subintervals
c     ppqad - integrates the pp-representation
c     bfqad - integrates the product of a function f and any spline
c             derivative in the b-representation
c     pfqad - integrates the product of a function f and any spline
c             derivative in the pp-representation
c     bvalu - evaluates the b-representation or a derivative
c     ppval - evaluates the pp-representation or a derivative
c     intrv - gets the largest index of the knot to the left of x
c     bsppp - converts from b- to pp-representation
c     bspvd - computes nonzero basis functions and derivatives at x
c     bspdr - sets up difference array for bspev
c     bspev - evaluates the b-representation and derivatives
c     bspvn - called by bspev, bspvd, bsppp and bintk for function and
c             derivative evaluations
c                        auxiliary routines
c
c       bsgq8,ppgq8,bnslv,bnfac,xermsg,dbsgq8,dppgq8,dbnslv,dbnfac
c
c                    machine dependent routines
c
c                      i1mach, r1mach, d1mach
c
c***references  1. d. e. amos, computation with splines and
c                 b-splines, report sand78-1968, sandia
c                 laboratories, march 1979.
c               2. d. e. amos, quadrature subroutines for splines and
c                 b-splines, report sand79-1825, sandia laboratories,
c                 december 1979.
c               3. carl de boor, a practical guide to splines, applied
c                 mathematics series 27, springer-verlag, new york,
c                 1978.
c               4. carl de boor, on calculating with b-splines, journal
c                 of approximation theory 6, (1972), pp. 50-62.
c               5. carl de boor, package for calculating with b-splines,
c                 siam journal on numerical analysis 14, 3 (june 1977),
c                 pp. 441-472.
c               6. r. j. hanson, constrained least squares curve fitting
c                 to discrete data using b-splines, a users guide,
c                 report sand78-1291, sandia laboratories, december
c                 1978.
c               7. f. n. fritsch and r. e. carlson, monotone piecewise
c                 cubic interpolation, siam journal on numerical ana-
c                 lysis 17, 2 (april 1980), pp. 238-246.
c***routines called  (none)
c***revision history  (yymmdd)
c   810223  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900723  purpose section revised.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  bspdoc
c***first executable statement  bspdoc
      return
      end
