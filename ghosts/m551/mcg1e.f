*deck %W%  %G%
      subroutine mcg1e(ni,nmi,loci,leni,mix,cm,r,g)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit real*8(a-h,o-z)
cc
cc
      dimension loci(2),leni(2),mix(2),cm(nmi,ni)
      dimension g(2)
      dimension r(nmi,ni)
c
      common /io/ inp,iout
c
c-----------------------------------------------------------------------
c
c --- description     this routine makes contributions from a block
c                     of 1-electron integrals to g(ij).
c                     the integrals (i j) are stored in array r.
c
c --- input
c
c     ni              no of i active orbitals.
c     nmi             no of i orbitals, including core and virtual.
c     loci(ni)        loci(i) starting position of the ith vector in
c                      the arrays mix and cm.
c     leni(ni)        length of the ith vector.
c     mix(--)         indices of vector components.
c     cm(--)          vector components.
c     r(nmi,ni)       transformed one-electron integrals.
c
c --- output
c
c     g(2)
c
c-----------------------------------------------------------------------
c
c      do 900 i = 1, ni
c 900 write (iout,9000) (r(j,i),j=1,nmi)
 9000 format(' *mcg1e r '//4(1x,f16.8))
c     write (iout,9001) (loci(i),i=1,ni)
 9001 format(' *mcg1e loc'//10(1x,i5))
c     write (iout,9002) (leni(i),i=1,ni)
 9002 format(' *mcg1e len'//10(1x,i5))
      ij = 1
c
      do 300 i = 1, ni
         mi1 = loci(i) + 1
         mi2 = mi1 + leni(i) - 1
cc
         do 400 j = 1, i
            mj1 = loci(j) + 1
            mj2 = mj1 + leni(j) - 1
            t = 0.d0
            if (mi2 .lt. mi1) go to 452
ccc
            do 450 m = mi1, mi2
               t = t + r(mix(m),j) * cm(mix(m),i)
c     write (iout,451) i, j, mix(m), r(mix(m),j), cm(mix(m),j), t
c 451 format(' *mcg1e i j mix r cm t ',3(1x,i5),2(1x,f16.8))
 450        continue
ccc
 452        if (j .ne. i) go to 455
            t = t + t
            go to 480
ccc
 455        if (mj2 .lt. mj1) go to 480
            do 460 m = mj1, mj2
               t = t + r(mix(m),i) * cm(mix(m),j)
c     write (iout,451) i, j, mix(m), r(mix(m),i), cm(mix(m),i), t
 460        continue
ccc
 480        g(ij) = g(ij) + t
            ij = ij + 1
 400     continue
cc
 300  continue
c
c     nij = ni * (ni + 1) / 2
c     write (iout,1002) (g(ij),ij=1,nij)
c1002 format(' *mcg1e g '//4(1x,f16.8))
      return
      end
