*deck @(#)mccg1e.f	5.1  11/6/94
      subroutine mccg1e(ni,nmi,cg,r,g)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mccg1e.f	5.1   11/6/94
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
      dimension g(2),r(nmi,ni)
      dimension cg(nmi,ni)
c
      common /io/ inp,iout
c
c-----------------------------------------------------------------------
c
c --- description     this routine makes contributions from a block
c                     of 1-electron integrals and the corresponding
c                     density matrix to c*g.
c
c --- input
c
c     ni              no of i active orbitals.
c     nmi             no of i orbitals, including core and virtual.
c     g(--)           density matrix elements.
c     r(nmi,ni)       transformed one-electron integrals.
c
c --- output
c
c     cg(nmi,ni)
c
c-----------------------------------------------------------------------
c     nij = ni * (ni + 1) / 2
c     write (iout,1002) (g(ij),ij=1,nij)
c1002 format(' *mccg1e g '//4(1x,f16.8))
c
      ij = 1
c
      do 300 i = 1, ni
cc
         do 400 j = 1, i
            t = g(ij)
            ij = ij + 1
ccc
            do 450 m = 1, nmi
               cg(m,i) = cg(m,i) + t * r(m,j)
               cg(m,j) = cg(m,j) + t * r(m,i)
 450        continue
ccc
 400     continue
cc
 300  continue
c
c     do 900 i = 1, ni
c 900 write (iout,9000) (cg(j,i),j=1,nmi)
c9000 format(' *mccg1e cg '//4(1x,f16.8))
c
      return
      end
