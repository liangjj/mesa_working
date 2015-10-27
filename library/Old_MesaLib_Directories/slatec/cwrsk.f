*deck cwrsk
      subroutine cwrsk (zr, fnu, kode, n, y, nz, cw, tol, elim, alim)
c***begin prologue  cwrsk
c***subsidiary
c***purpose  subsidiary to cbesi and cbesk
c***library   slatec
c***type      all (cwrsk-a, zwrsk-a)
c***author  amos, d. e., (snl)
c***description
c
c     cwrsk computes the i bessel function for re(z).ge.0.0 by
c     normalizing the i function ratios from crati by the wronskian
c
c***see also  cbesi, cbesk
c***routines called  cbknu, crati, r1mach
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  cwrsk
      complex cinu, cscl, ct, cw, c1, c2, rct, st, y, zr
      real act, acw, alim, ascle, elim, fnu, s1, s2, tol, yy, r1mach
      integer i, kode, n, nw, nz
      dimension y(n), cw(2)
c***first executable statement  cwrsk
c-----------------------------------------------------------------------
c     i(fnu+i-1,z) by backward recurrence for ratios
c     y(i)=i(fnu+i,z)/i(fnu+i-1,z) from crati normalized by the
c     wronskian with k(fnu,z) and k(fnu+1,z) from cbknu.
c-----------------------------------------------------------------------
      nz = 0
      call cbknu(zr, fnu, kode, 2, cw, nw, tol, elim, alim)
      if (nw.ne.0) go to 50
      call crati(zr, fnu, n, y, tol)
c-----------------------------------------------------------------------
c     recur forward on i(fnu+1,z) = r(fnu,z)*i(fnu,z),
c     r(fnu+j-1,z)=y(j),  j=1,...,n
c-----------------------------------------------------------------------
      cinu = cmplx(1.0e0,0.0e0)
      if (kode.eq.1) go to 10
      yy = aimag(zr)
      s1 = cos(yy)
      s2 = sin(yy)
      cinu = cmplx(s1,s2)
   10 continue
c-----------------------------------------------------------------------
c     on low exponent machines the k functions can be close to both
c     the under and overflow limits and the normalization must be
c     scaled to prevent over or underflow. cuoik has determined that
c     the result is on scale.
c-----------------------------------------------------------------------
      acw = abs(cw(2))
      ascle = 1.0e+3*r1mach(1)/tol
      cscl = cmplx(1.0e0,0.0e0)
      if (acw.gt.ascle) go to 20
      cscl = cmplx(1.0e0/tol,0.0e0)
      go to 30
   20 continue
      ascle = 1.0e0/ascle
      if (acw.lt.ascle) go to 30
      cscl = cmplx(tol,0.0e0)
   30 continue
      c1 = cw(1)*cscl
      c2 = cw(2)*cscl
      st = y(1)
c-----------------------------------------------------------------------
c     cinu=cinu*(conjg(ct)/abs(ct))*(1.0e0/abs(ct) prevents
c     under- or overflow prematurely by squaring abs(ct)
c-----------------------------------------------------------------------
      ct = zr*(c2+st*c1)
      act = abs(ct)
      rct = cmplx(1.0e0/act,0.0e0)
      ct = conjg(ct)*rct
      cinu = cinu*rct*ct
      y(1) = cinu*cscl
      if (n.eq.1) return
      do 40 i=2,n
        cinu = st*cinu
        st = y(i)
        y(i) = cinu*cscl
   40 continue
      return
   50 continue
      nz = -1
      if(nw.eq.(-2)) nz=-2
      return
      end
