*deck %W%  %G%
      subroutine tag1e(nao,nob,ta,r,g,ndf,sg)
cc
cc
      implicit real*8(a-h,o-z)
      dimension ta(*)
      dimension g(*),sg(*)
      dimension r(*)
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
c     nao              no of i active orbitals.
c     nob             no of i orbitals, including core and virtual.
c     loci(ni)        loci(i) starting position of the ith vector in
c                      the arrays mix and cm.
c     leni(ni)        length of the ith vector.
c     mix(--)         indices of vector components.
c     cm(--)          vector components.
c     r(nob,nao)       transformed one-electron integrals.
c
c --- output
c
c     g(2)
c
c-----------------------------------------------------------------------
c
c
      mk=nao*ndf
      nnao=nao*(nao+1)/2
      naonao=nao*nao
c
c..bhl
c      write(iout,*)' tag1e:  ta(nob,nao)  mo ndf=1 '
c      call matout(ta,nob,nao,nob,nao,iout)
c      write(iout,*)'  tag1e: r(nob,nao)  mo  '
c      call matout(r,nob,nao,nob,nao,iout)
c..bhl
c
      call ebtc(sg,r,ta,nao,nob,mk)
c
      ix=1
      jx=1
c..bhl
c      write(iout,*)' tag1e: before fold ta_g1 ndf=1 '
c      call matout(sg,nao,nao,nao,iout)
c..bhl
      do 10 i=1,ndf
         call blfold(g(ix),sg(jx),nao)
         ix=ix+nnao
         jx=jx+naonao
  10  continue
c..bhl
c      write(iout,*)'  ta_g1 ndf=1 '
c      write(iout,11)(g(ii),ii=1,nnao)
c  11  format(5(1x,f12.8))
c..bhl
c
c
      return
      end
