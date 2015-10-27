*deck bkisr
      subroutine bkisr (x, n, sum, ierr)
c***begin prologue  bkisr
c***subsidiary
c***purpose  subsidiary to bskin
c***library   slatec
c***type      single precision (bkisr-s, dbkisr-d)
c***author  amos, d. e., (snla)
c***description
c
c     bkisr computes repeated integrals of the k0 bessel function
c     by the series for n=0,1, and 2.
c
c***see also  bskin
c***routines called  psixn, r1mach
c***revision history  (yymmdd)
c   820601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  bkisr
      integer i, ierr, k, kk, kkn, k1, n, np
      real ak, atol, bk, c, fk, fn, hx, hxs, pol, pr, sum, tkp, tol,
     * trm, x, xln
      real psixn, r1mach
      dimension c(2)
      save c
c
      data c(1), c(2) /1.57079632679489662e+00,1.0e0/
c***first executable statement  bkisr
      ierr=0
      tol = max(r1mach(4),1.0e-18)
      if (x.lt.tol) go to 50
      pr = 1.0e0
      pol = 0.0e0
      if (n.eq.0) go to 20
      do 10 i=1,n
        pol = -pol*x + c(i)
        pr = pr*x/i
   10 continue
   20 continue
      hx = x*0.5e0
      hxs = hx*hx
      xln = log(hx)
      np = n + 1
      tkp = 3.0e0
      fk = 2.0e0
      fn = n
      bk = 4.0e0
      ak = 2.0e0/((fn+1.0e0)*(fn+2.0e0))
      sum = ak*(psixn(n+3)-psixn(3)+psixn(2)-xln)
      atol = sum*tol*0.75e0
      do 30 k=2,20
        ak = ak*(hxs/bk)*((tkp+1.0e0)/(tkp+fn+1.0e0))*(tkp/(tkp+fn))
        k1 = k + 1
        kk = k1 + k
        kkn = kk + n
        trm = (psixn(k1)+psixn(kkn)-psixn(kk)-xln)*ak
        sum = sum + trm
        if (abs(trm).le.atol) go to 40
        tkp = tkp + 2.0e0
        bk = bk + tkp
        fk = fk + 1.0e0
   30 continue
      go to 80
   40 continue
      sum = (sum*hxs+psixn(np)-xln)*pr
      if (n.eq.1) sum = -sum
      sum = pol + sum
      return
c-----------------------------------------------------------------------
c     small x case, x.lt.word tolerance
c-----------------------------------------------------------------------
   50 continue
      if (n.gt.0) go to 60
      hx = x*0.5e0
      sum = psixn(1) - log(hx)
      return
   60 continue
      sum = c(n)
      return
   80 continue
      ierr=2
      return
      end
