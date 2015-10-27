*deck dbkisr
      subroutine dbkisr (x, n, sum, ierr)
c***begin prologue  dbkisr
c***subsidiary
c***purpose  subsidiary to dbskin
c***library   slatec
c***type      double precision (bkisr-s, dbkisr-d)
c***author  amos, d. e., (snla)
c***description
c
c     dbkisr computes repeated integrals of the k0 bessel function
c     by the series for n=0,1, and 2.
c
c***see also  dbskin
c***routines called  d1mach, dpsixn
c***revision history  (yymmdd)
c   820601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  dbkisr
      integer i, ierr, k, kk, kkn, k1, n, np
      double precision ak, atol, bk, c, fk, fn, hx, hxs, pol, pr, sum,
     * tkp, tol, trm, x, xln
      double precision dpsixn, d1mach
      dimension c(2)
      save c
c
      data c(1), c(2) /1.57079632679489662d+00,1.0d0/
c***first executable statement  dbkisr
      ierr=0
      tol = max(d1mach(4),1.0d-18)
      if (x.lt.tol) go to 50
      pr = 1.0d0
      pol = 0.0d0
      if (n.eq.0) go to 20
      do 10 i=1,n
        pol = -pol*x + c(i)
        pr = pr*x/i
   10 continue
   20 continue
      hx = x*0.5d0
      hxs = hx*hx
      xln = log(hx)
      np = n + 1
      tkp = 3.0d0
      fk = 2.0d0
      fn = n
      bk = 4.0d0
      ak = 2.0d0/((fn+1.0d0)*(fn+2.0d0))
      sum = ak*(dpsixn(n+3)-dpsixn(3)+dpsixn(2)-xln)
      atol = sum*tol*0.75d0
      do 30 k=2,20
        ak = ak*(hxs/bk)*((tkp+1.0d0)/(tkp+fn+1.0d0))*(tkp/(tkp+fn))
        k1 = k + 1
        kk = k1 + k
        kkn = kk + n
        trm = (dpsixn(k1)+dpsixn(kkn)-dpsixn(kk)-xln)*ak
        sum = sum + trm
        if (abs(trm).le.atol) go to 40
        tkp = tkp + 2.0d0
        bk = bk + tkp
        fk = fk + 1.0d0
   30 continue
      go to 80
   40 continue
      sum = (sum*hxs+dpsixn(np)-xln)*pr
      if (n.eq.1) sum = -sum
      sum = pol + sum
      return
c-----------------------------------------------------------------------
c     small x case, x.lt.word tolerance
c-----------------------------------------------------------------------
   50 continue
      if (n.gt.0) go to 60
      hx = x*0.5d0
      sum = dpsixn(1) - log(hx)
      return
   60 continue
      sum = c(n)
      return
   80 continue
      ierr=2
      return
      end
