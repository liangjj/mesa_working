*deck passb4
      subroutine passb4 (ido, l1, cc, ch, wa1, wa2, wa3)
c***begin prologue  passb4
c***subsidiary
c***purpose  calculate the fast fourier transform of subvectors of
c            length four.
c***library   slatec (fftpack)
c***type      single precision (passb4-s)
c***author  swarztrauber, p. n., (ncar)
c***routines called  (none)
c***revision history  (yymmdd)
c   790601  date written
c   830401  modified to use slatec library source file format.
c   860115  modified by ron boisvert to adhere to fortran 77 by
c           changing dummy array size declarations (1) to (*).
c   881128  modified by dick valent to meet prologue standards.
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  passb4
      dimension cc(ido,4,*), ch(ido,l1,4), wa1(*), wa2(*), wa3(*)
c***first executable statement  passb4
      if (ido .ne. 2) go to 102
      do 101 k=1,l1
         ti1 = cc(2,1,k)-cc(2,3,k)
         ti2 = cc(2,1,k)+cc(2,3,k)
         tr4 = cc(2,4,k)-cc(2,2,k)
         ti3 = cc(2,2,k)+cc(2,4,k)
         tr1 = cc(1,1,k)-cc(1,3,k)
         tr2 = cc(1,1,k)+cc(1,3,k)
         ti4 = cc(1,2,k)-cc(1,4,k)
         tr3 = cc(1,2,k)+cc(1,4,k)
         ch(1,k,1) = tr2+tr3
         ch(1,k,3) = tr2-tr3
         ch(2,k,1) = ti2+ti3
         ch(2,k,3) = ti2-ti3
         ch(1,k,2) = tr1+tr4
         ch(1,k,4) = tr1-tr4
         ch(2,k,2) = ti1+ti4
         ch(2,k,4) = ti1-ti4
  101 continue
      return
  102 if(ido/2.lt.l1) go to 105
      do 104 k=1,l1
cdir$ ivdep
         do 103 i=2,ido,2
            ti1 = cc(i,1,k)-cc(i,3,k)
            ti2 = cc(i,1,k)+cc(i,3,k)
            ti3 = cc(i,2,k)+cc(i,4,k)
            tr4 = cc(i,4,k)-cc(i,2,k)
            tr1 = cc(i-1,1,k)-cc(i-1,3,k)
            tr2 = cc(i-1,1,k)+cc(i-1,3,k)
            ti4 = cc(i-1,2,k)-cc(i-1,4,k)
            tr3 = cc(i-1,2,k)+cc(i-1,4,k)
            ch(i-1,k,1) = tr2+tr3
            cr3 = tr2-tr3
            ch(i,k,1) = ti2+ti3
            ci3 = ti2-ti3
            cr2 = tr1+tr4
            cr4 = tr1-tr4
            ci2 = ti1+ti4
            ci4 = ti1-ti4
            ch(i-1,k,2) = wa1(i-1)*cr2-wa1(i)*ci2
            ch(i,k,2) = wa1(i-1)*ci2+wa1(i)*cr2
            ch(i-1,k,3) = wa2(i-1)*cr3-wa2(i)*ci3
            ch(i,k,3) = wa2(i-1)*ci3+wa2(i)*cr3
            ch(i-1,k,4) = wa3(i-1)*cr4-wa3(i)*ci4
            ch(i,k,4) = wa3(i-1)*ci4+wa3(i)*cr4
  103    continue
  104 continue
      return
  105 do 107 i=2,ido,2
cdir$ ivdep
         do 106 k=1,l1
            ti1 = cc(i,1,k)-cc(i,3,k)
            ti2 = cc(i,1,k)+cc(i,3,k)
            ti3 = cc(i,2,k)+cc(i,4,k)
            tr4 = cc(i,4,k)-cc(i,2,k)
            tr1 = cc(i-1,1,k)-cc(i-1,3,k)
            tr2 = cc(i-1,1,k)+cc(i-1,3,k)
            ti4 = cc(i-1,2,k)-cc(i-1,4,k)
            tr3 = cc(i-1,2,k)+cc(i-1,4,k)
            ch(i-1,k,1) = tr2+tr3
            cr3 = tr2-tr3
            ch(i,k,1) = ti2+ti3
            ci3 = ti2-ti3
            cr2 = tr1+tr4
            cr4 = tr1-tr4
            ci2 = ti1+ti4
            ci4 = ti1-ti4
            ch(i-1,k,2) = wa1(i-1)*cr2-wa1(i)*ci2
            ch(i,k,2) = wa1(i-1)*ci2+wa1(i)*cr2
            ch(i-1,k,3) = wa2(i-1)*cr3-wa2(i)*ci3
            ch(i,k,3) = wa2(i-1)*ci3+wa2(i)*cr3
            ch(i-1,k,4) = wa3(i-1)*cr4-wa3(i)*ci4
            ch(i,k,4) = wa3(i-1)*ci4+wa3(i)*cr4
  106    continue
  107 continue
      return
      end
