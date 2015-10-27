*deck @(#)mccgpk.f	5.1  11/6/94
      subroutine mccgpk(nobj,nobi,cgx,lok,len,mix,cg,l)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mccgpk.f	5.1   11/6/94
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
      dimension lok(2),len(2),mix(2),cgx(nobj,nobi)
      dimension cg(2)
c
      common /io/ inp,iout
c
c
c-----------------------------------------------------------------------
c
c --- description     this program makes contributions from a block
c                     of 1-symmetry integrals to g(ijkl).  the integrals
c                     (ij kl) are stored in canonical order.
c
c --- input
c
c     nobj            no of orbitals.
c     nobi            no of occupied orbitals.
c     cgx(nobj,nobi)  c*g expanded form.
c     lok(nobi)       loc(k) starting postion of the kth vector in
c                      the arrays mix and cm.
c     len(nobi)       length of the kth vector.
c     mix(--)         indices of vector components.
c
c --- output
c
c     cg(2)
c
c-----------------------------------------------------------------------
c     write (iout,9876) (lok(i),i=1,nobi)
c9876 format(' *mccgpk lok '//10(1x,i4))
c     write (iout,9875) (len(i),i=1,nobi)
c9875 format(' *mccgpk len '//10(1x,i4))
c
      do 100 i = 1, nobi
         j1 = i + 1
         if (j1 .gt. nobi) go to 250
cc
         do 200 j = j1, nobi
 200     cgx(j,i) = cgx(i,j) - cgx(j,i)
cc
 250     j1 = nobi + 1
         if (j1 .gt. nobj) go to 100
cc
         do 300 j = j1, nobj
 300     cgx(j,i) = - cgx(j,i)
cc
 100  continue
c
      l = 0
      do 500 i = 1, nobi
         m1 = lok(i) + 1
         m2 = m1 + len(i) - 1
         if (m2 .lt. m1) go to 500
cc
         do 600 m = m1, m2
c     write (iout,601) i, m, mix(m), l
c 601 format(' *mccgpk i m mix(m) l',4(1x,i6))
c
            if (i .ge. mix(m)) go to 600
            l = l + 1
            cg(l) = cg(l) - cgx(mix(m),i)
 600     continue
cc
 500  continue
      return
      end
