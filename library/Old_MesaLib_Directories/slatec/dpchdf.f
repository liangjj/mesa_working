*deck dpchdf
      double precision function dpchdf (k, x, s, ierr)
c***begin prologue  dpchdf
c***subsidiary
c***purpose  computes divided differences for dpchce and dpchsp
c***library   slatec (pchip)
c***type      double precision (pchdf-s, dpchdf-d)
c***author  fritsch, f. n., (llnl)
c***description
c
c          dpchdf:   dpchip finite difference formula
c
c     uses a divided difference formulation to compute a k-point approx-
c     imation to the derivative at x(k) based on the data in x and s.
c
c     called by  dpchce  and  dpchsp  to compute 3- and 4-point boundary
c     derivative approximations.
c
c ----------------------------------------------------------------------
c
c     on input:
c        k      is the order of the desired derivative approximation.
c               k must be at least 3 (error return if not).
c        x      contains the k values of the independent variable.
c               x need not be ordered, but the values **must** be
c               distinct.  (not checked here.)
c        s      contains the associated slope values:
c                  s(i) = (f(i+1)-f(i))/(x(i+1)-x(i)), i=1(1)k-1.
c               (note that s need only be of length k-1.)
c
c     on return:
c        s      will be destroyed.
c        ierr   will be set to -1 if k.lt.2 .
c        dpchdf  will be set to the desired derivative approximation if
c               ierr=0 or to zero if ierr=-1.
c
c ----------------------------------------------------------------------
c
c***see also  dpchce, dpchsp
c***references  carl de boor, a practical guide to splines, springer-
c                 verlag, new york, 1978, pp. 10-16.
c***routines called  xermsg
c***revision history  (yymmdd)
c   820503  date written
c   820805  converted to slatec library version.
c   870707  corrected xerror calls for d.p. name(s).
c   870813  minor cosmetic changes.
c   890206  corrected xerror calls.
c   890411  added save statements (vers. 3.2).
c   890411  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900328  added type section.  (wrb)
c   910408  updated author and date written sections in prologue.  (wrb)
c   920429  revised format and order of references.  (wrb,fnf)
c   930503  improved purpose.  (fnf)
c***end prologue  dpchdf
c
c**end
c
c  declare arguments.
c
      integer  k, ierr
      double precision  x(k), s(k)
c
c  declare local variables.
c
      integer  i, j
      double precision  value, zero
      save zero
      data  zero /0.d0/
c
c  check for legal value of k.
c
c***first executable statement  dpchdf
      if (k .lt. 3)  go to 5001
c
c  compute coefficients of interpolating polynomial.
c
      do 10  j = 2, k-1
         do 9  i = 1, k-j
            s(i) = (s(i+1)-s(i))/(x(i+j)-x(i))
    9    continue
   10 continue
c
c  evaluate derivative at x(k).
c
      value = s(1)
      do 20  i = 2, k-1
         value = s(i) + value*(x(k)-x(i))
   20 continue
c
c  normal return.
c
      ierr = 0
      dpchdf = value
      return
c
c  error return.
c
 5001 continue
c     k.lt.3 return.
      ierr = -1
      call xermsg ('slatec', 'dpchdf', 'k less than three', ierr, 1)
      dpchdf = zero
      return
c------------- last line of dpchdf follows -----------------------------
      end
