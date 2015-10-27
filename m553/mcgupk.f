*deck @(#)mcgupk.f	5.1  11/6/94
      subroutine mcgupk(nobj,nobi,ncob,cg,lok,len,mix,cgx,l)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcgupk.f	5.1   11/6/94
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
      dimension lok(2),len(2),mix(2),cg(2)
      dimension cgx(nobj,nobi)
c
      common /io/ inp,iout
c
c
c-----------------------------------------------------------------------
c
c --- description     this routine makes contributions from a block
c                     of 1-symmetry integrals to g(ijkl).  the integrals
c                     (ij kl) are stored in canonical order.
c
c --- input
c
c     nobj            no of orbitals.
c     nobi            no of occupied orbitals.
c     cgx(nobj,nobi)  c*g expanded form.
c     loc(nobi)       loc(k) starting postion of the kth vector in
c                      the arrays mix and cm.
c     len(nobi)       length of the kth vector.
c     mix(--)         indices of vector components.
c
c --- output
c
c     cg(2)
c
c-----------------------------------------------------------------------
c     write (iout,9876) (loc(i),i=1,nobi)
c9876 format(' *mcgupk loc '//10(1x,i4))
c     write (iout,9875) (len(i),i=1,nobi)
c9875 format(' *mcgupk len '//10(1x,i4))
 
      do 50 i = 1, nobi
         do 55 j = 1, nobj
 55      cgx(j,i) = 0.d0
 50   continue
c
      l = 0
      do 500 i = 1, nobi
         m1 = lok(i) + 1
         m2 = m1 + len(i) - 1
         if (m2 .lt. m1) go to 500
cc
         do 600 m = m1, m2
c     write (iout,601) i, m, mix(m)
c 601 format(' *mcgupk i m mix(m)',3(1x,i6))
c
 605        if (i .ge. mix(m)) go to 600
            l = l + 1
            cgx(mix(m),i) = cg(l)
 600     continue
cc
 500  continue
c
      do 100 i = 1, nobi
         j1 = i - 1
         if (j1 .lt. 1) go to 100
cc
         do 200 j = 1, j1
 200     cgx(j,i) = - cgx(i,j)
cc
 100  continue
c
c     do 800 i = 1, nobi
c     write (iout,801) (cgx(j,i),j=1,nobj)
c 800 continue
c 801 format(' *mcgupk cgx '//4(1x,f16.8))
c     write (iout,802) (cg(ll),ll=1,l)
c 802 format(' *mcgupk cg '//4(1x,f16.8))
      return
      end
