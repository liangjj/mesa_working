*deck fzero
c***begin prologue  fzero
c***purpose  search for a zero of a function f(x) in a given interval
c            (b,c).  it is designed primarily for problems where f(b)
c            and f(c) have opposite signs.
c***library   slatec
c***category  f1b
c***type      double precision (fzero-s, dfzero-d)
c***keywords  bisection, nonlinear, roots, zeros
c***author  shampine, l. f., (snla)
c           watts, h. a., (snla)
c***description
c
c     fzero searches for a zero of a double precision function f(x)
c     between the given double precision values b and c until the width
c     of the interval (b,c) has collapsed to within a tolerance
c     specified by the stopping criterion,
c        abs(b-c) .le. 2.*(rw*abs(b)+ae).
c     the method used is an efficient combination of bisection and the
c     secant rule and is due to t. j. dekker.
c
c     description of arguments
c
c   f     :ext   - name of the double precision external function.  this
c                  name must be in an external statement in the calling
c                  program.  f must be a function of one double
c                  precision argument.
c
c   b     :inout - one end of the double precision interval (b,c).  the
c                  value returned for b usually is the better
c                  approximation to a zero of f.
c
c   c     :inout - the other end of the double precision interval (b,c)
c
c   r     :in    - a (better) double precision guess of a zero of f
c                  which could help in speeding up convergence.  if f(b)
c                  and f(r) have opposite signs, a root will be found in
c                  the interval (b,r);  if not, but f(r) and f(c) have
c                  opposite signs, a root will be found in the interval
c                  (r,c);  otherwise, the interval (b,c) will be
c                  searched for a possible root.  when no better guess
c                  is known, it is recommended that r be set to b or c,
c                  since if r is not interior to the interval (b,c), it
c                  will be ignored.
c
c   re    :in    - relative error used for rw in the stopping criterion.
c                  if the requested re is less than machine precision,
c                  then rw is set to approximately machine precision.
c
c   ae    :in    - absolute error used in the stopping criterion.  if
c                  the given interval (b,c) contains the origin, then a
c                  nonzero value should be chosen for ae.
c
c   iflag :out   - a status code.  user must check iflag after each
c                  call.  control returns to the user from dfzero in all
c                  cases.
c
c                1  b is within the requested tolerance of a zero.
c                   the interval (b,c) collapsed to the requested
c                   tolerance, the function changes sign in (b,c), and
c                   f(x) decreased in magnitude as (b,c) collapsed.
c
c                2  f(b) = 0.  however, the interval (b,c) may not have
c                   collapsed to the requested tolerance.
c
c                3  b may be near a singular point of f(x).
c                   the interval (b,c) collapsed to the requested tol-
c                   erance and the function changes sign in (b,c), but
c                   f(x) increased in magnitude as (b,c) collapsed, i.e.
c                     abs(f(b out)) .gt. max(abs(f(b in)),abs(f(c in)))
c
c                4  no change in sign of f(x) was found although the
c                   interval (b,c) collapsed to the requested tolerance.
c                   the user must examine this case and decide whether
c                   b is near a local minimum of f(x), or b is near a
c                   zero of even multiplicity, or neither of these.
c
c                5  too many (.gt. maxit) function evaluations used.
c
c***references  l. f. shampine and h. a. watts, fzero, a root-solving
c                 code, report sc-tm-70-631, sandia laboratories,
c                 september 1970.
c               t. j. dekker, finding a zero by means of successive
c                 linear interpolation, constructive aspects of the
c                 fundamental theorem of algebra, edited by b. dejon
c                 and p. henrici, wiley-interscience, 1969.
c***routines called  d1mach
c***revision history  (yymmdd)
c   700901  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dfzero
      subroutine fzero ( b, c, r, re, ae, maxit,iflag,
     1                  band,point,v,ipiv,work,stp,n,order,bw,
     2                  bc,type)
      real*8 a,acbs,acmb,ae,aw,b,c,cmb,r1mach,er,
     +         fa,fb,fc,fx,fz,p,q,r,re,rw,t,tol,z
      integer ic,iflag,kount
      real*8 band, point, v, work, stp
      character*(*) type, bc
      dimension band(n,-bw:bw), point(n), v(n) ,work(n,3), ipiv(n,2)
c
c***first executable statement  dfzero
c
c   er is two times the computer unit roundoff value which is defined
c   here by the function d1mach.
c
      er = 2.0d0 * r1mach(4)
c
c   initialize.
c
      z = r
      if (r .le. min(b,c)  .or.  r .ge. max(b,c)) z = c
      rw = max(re,er)
      aw = max(ae,0.d0)
      ic = 0
      t = z
      call condit(t,fz,band,point,v,ipiv,work,stp,n,order,bw,bc,type)
      fc = fz
      t = b
      call condit(t,fb,band,point,v,ipiv,work,stp,n,order,bw,bc,type)
      kount = 2
      if (sign(1.0d0,fz) .eq. sign(1.0d0,fb)) go to 1
      c = z
      go to 2
    1 if (z .eq. c) go to 2
      t = c
      call condit(t,fc,band,point,v,ipiv,work,stp,n,order,bw,bc,type)
      kount = 3
      if (sign(1.0d0,fz) .eq. sign(1.0d0,fc)) go to 2
      b = z
      fb = fz
    2 a = c
      fa = fc
      acbs = abs(b-c)
      fx = max(abs(fb),abs(fc))
c
    3 if (abs(fc) .ge. abs(fb)) go to 4
c
c   perform interchange.
c
      a = b
      fa = fb
      b = c
      fb = fc
      c = a
      fc = fa
c
    4 cmb = 0.5d0*(c-b)
      acmb = abs(cmb)
      tol = rw*abs(b) + aw
c
c   test stopping criterion and function count.
c
      if (acmb .le. tol) go to 10
      if (fb .eq. 0.d0) go to 11
      if (kount .ge. maxit) go to 14
c
c   calculate new iterate implicitly as b+p/q, where we arrange
c   p .ge. 0.  the implicit form is used to prevent overflow.
c
      p = (b-a)*fb
      q = fa - fb
      if (p .ge. 0.d0) go to 5
      p = -p
      q = -q
c
c   update a and check for satisfactory reduction in the size of the
c   bracketing interval.  if not, perform bisection.
c
    5 a = b
      fa = fb
      ic = ic + 1
      if (ic .lt. 4) go to 6
      if (8.0d0*acmb .ge. acbs) go to 8
      ic = 0
      acbs = acmb
c
c   test for too small a change.
c
    6 if (p .gt. abs(q)*tol) go to 7
c
c   increment by tolerance.
c
      b = b + sign(tol,cmb)
      go to 9
c
c   root ought to be between b and (c+b)/2.
c
    7 if (p .ge. cmb*q) go to 8
c
c   use secant rule.
c
      b = b + p/q
      go to 9
c
c   use bisection (c+b)/2.
c
    8 b = b + cmb
c
c   have completed computation for new iterate b.
c
    9 t = b
      call condit(t,fb,band,point,v,ipiv,work,stp,n,order,bw,bc,type)
      kount = kount + 1
c
c   decide whether next step is interpolation or extrapolation.
c
      if (sign(1.0d0,fb) .ne. sign(1.0d0,fc)) go to 3
      c = a
      fc = fa
      go to 3
c
c   finished.  process results for proper setting of iflag.
c
   10 if (sign(1.0d0,fb) .eq. sign(1.0d0,fc)) go to 13
      if (abs(fb) .gt. fx) go to 12
      iflag = 1
      return
   11 iflag = 2
      return
   12 iflag = 3
      return
   13 iflag = 4
      return
   14 iflag = 5
      return
      end
