*deck intyd
      subroutine intyd (t, k, yh, nyh, dky, iflag)
c***begin prologue  intyd
c***subsidiary
c***purpose  subsidiary to debdf
c***library   slatec
c***type      single precision (intyd-s, dintyd-d)
c***author  watts, h. a., (snla)
c***description
c
c   intyd approximates the solution and derivatives at t by polynomial
c   interpolation. must be used in conjunction with the integrator
c   package debdf.
c ----------------------------------------------------------------------
c intyd computes interpolated values of the k-th derivative of the
c dependent variable vector y, and stores it in dky.
c this routine is called by debdf with k = 0,1 and t = tout, but may
c also be called by the user for any k up to the current order.
c (see detailed instructions in lsode usage documentation.)
c ----------------------------------------------------------------------
c the computed values in dky are gotten by interpolation using the
c nordsieck history array yh.  this array corresponds uniquely to a
c vector-valued polynomial of degree nqcur or less, and dky is set
c to the k-th derivative of this polynomial at t.
c the formula for dky is..
c              q
c  dky(i)  =  sum  c(j,k) * (t - tn)**(j-k) * h**(-j) * yh(i,j+1)
c             j=k
c where  c(j,k) = j*(j-1)*...*(j-k+1), q = nqcur, tn = tcur, h = hcur.
c the quantities  nq = nqcur, l = nq+1, n = neq, tn, and h are
c communicated by common.  the above sum is done in reverse order.
c iflag is returned negative if either k or t is out of bounds.
c ----------------------------------------------------------------------
c
c***see also  debdf
c***routines called  (none)
c***common blocks    debdf1
c***revision history  (yymmdd)
c   800901  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  intyd
c
clll. optimize
      integer k, nyh, iflag, i, ic, ier, iownd, iowns, j, jb, jb2,
     1   jj, jj1, jp1, jstart, kflag, l, maxord, meth, miter, n, nfe,
     2   nje, nq, nqu, nst
      real t, yh, dky,
     1   rownd, rowns, el0, h, hmin, hmxi, hu, tn, uround,
     2   c, r, s, tp
      dimension yh(nyh,*), dky(*)
      common /debdf1/ rownd, rowns(210),
     1   el0, h, hmin, hmxi, hu, tn, uround, iownd(14), iowns(6),
     2   ier, jstart, kflag, l, meth, miter, maxord, n, nq, nst, nfe,
     3   nje, nqu
c
c***first executable statement  intyd
      iflag = 0
      if (k .lt. 0 .or. k .gt. nq) go to 80
      tp = tn - hu*(1.0e0 + 100.0e0*uround)
      if ((t-tp)*(t-tn) .gt. 0.0e0) go to 90
c
      s = (t - tn)/h
      ic = 1
      if (k .eq. 0) go to 15
      jj1 = l - k
      do 10 jj = jj1,nq
 10     ic = ic*jj
 15   c = ic
      do 20 i = 1,n
 20     dky(i) = c*yh(i,l)
      if (k .eq. nq) go to 55
      jb2 = nq - k
      do 50 jb = 1,jb2
        j = nq - jb
        jp1 = j + 1
        ic = 1
        if (k .eq. 0) go to 35
        jj1 = jp1 - k
        do 30 jj = jj1,j
 30       ic = ic*jj
 35     c = ic
        do 40 i = 1,n
 40       dky(i) = c*yh(i,jp1) + s*dky(i)
 50     continue
      if (k .eq. 0) return
 55   r = h**(-k)
      do 60 i = 1,n
 60     dky(i) = r*dky(i)
      return
c
 80   iflag = -1
      return
 90   iflag = -2
      return
c----------------------- end of subroutine intyd -----------------------
      end
