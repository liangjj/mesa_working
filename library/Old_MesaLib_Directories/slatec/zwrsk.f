*deck zwrsk
      subroutine zwrsk (zrr, zri, fnu, kode, n, yr, yi, nz, cwr, cwi,
     +   tol, elim, alim)
c***begin prologue  zwrsk
c***subsidiary
c***purpose  subsidiary to zbesi and zbesk
c***library   slatec
c***type      all (cwrsk-a, zwrsk-a)
c***author  amos, d. e., (snl)
c***description
c
c     zwrsk computes the i bessel function for re(z).ge.0.0 by
c     normalizing the i function ratios from zrati by the wronskian
c
c***see also  zbesi, zbesk
c***routines called  d1mach, zabs, zbknu, zrati
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  zwrsk
c     complex cinu,cscl,ct,cw,c1,c2,rct,st,y,zr
      double precision act, acw, alim, ascle, cinui, cinur, csclr, cti,
     * ctr, cwi, cwr, c1i, c1r, c2i, c2r, elim, fnu, pti, ptr, ract,
     * sti, str, tol, yi, yr, zri, zrr, zabs, d1mach
      integer i, kode, n, nw, nz
      dimension yr(n), yi(n), cwr(2), cwi(2)
      external zabs
c***first executable statement  zwrsk
c-----------------------------------------------------------------------
c     i(fnu+i-1,z) by backward recurrence for ratios
c     y(i)=i(fnu+i,z)/i(fnu+i-1,z) from crati normalized by the
c     wronskian with k(fnu,z) and k(fnu+1,z) from cbknu.
c-----------------------------------------------------------------------
c
      nz = 0
      call zbknu(zrr, zri, fnu, kode, 2, cwr, cwi, nw, tol, elim, alim)
      if (nw.ne.0) go to 50
      call zrati(zrr, zri, fnu, n, yr, yi, tol)
c-----------------------------------------------------------------------
c     recur forward on i(fnu+1,z) = r(fnu,z)*i(fnu,z),
c     r(fnu+j-1,z)=y(j),  j=1,...,n
c-----------------------------------------------------------------------
      cinur = 1.0d0
      cinui = 0.0d0
      if (kode.eq.1) go to 10
      cinur = cos(zri)
      cinui = sin(zri)
   10 continue
c-----------------------------------------------------------------------
c     on low exponent machines the k functions can be close to both
c     the under and overflow limits and the normalization must be
c     scaled to prevent over or underflow. cuoik has determined that
c     the result is on scale.
c-----------------------------------------------------------------------
      acw = zabs(cwr(2),cwi(2))
      ascle = 1.0d+3*d1mach(1)/tol
      csclr = 1.0d0
      if (acw.gt.ascle) go to 20
      csclr = 1.0d0/tol
      go to 30
   20 continue
      ascle = 1.0d0/ascle
      if (acw.lt.ascle) go to 30
      csclr = tol
   30 continue
      c1r = cwr(1)*csclr
      c1i = cwi(1)*csclr
      c2r = cwr(2)*csclr
      c2i = cwi(2)*csclr
      str = yr(1)
      sti = yi(1)
c-----------------------------------------------------------------------
c     cinu=cinu*(conjg(ct)/abs(ct))*(1.0d0/abs(ct) prevents
c     under- or overflow prematurely by squaring abs(ct)
c-----------------------------------------------------------------------
      ptr = str*c1r - sti*c1i
      pti = str*c1i + sti*c1r
      ptr = ptr + c2r
      pti = pti + c2i
      ctr = zrr*ptr - zri*pti
      cti = zrr*pti + zri*ptr
      act = zabs(ctr,cti)
      ract = 1.0d0/act
      ctr = ctr*ract
      cti = -cti*ract
      ptr = cinur*ract
      pti = cinui*ract
      cinur = ptr*ctr - pti*cti
      cinui = ptr*cti + pti*ctr
      yr(1) = cinur*csclr
      yi(1) = cinui*csclr
      if (n.eq.1) return
      do 40 i=2,n
        ptr = str*cinur - sti*cinui
        cinui = str*cinui + sti*cinur
        cinur = ptr
        str = yr(i)
        sti = yi(i)
        yr(i) = cinur*csclr
        yi(i) = cinui*csclr
   40 continue
      return
   50 continue
      nz = -1
      if(nw.eq.(-2)) nz=-2
      return
      end
