*deck @(#)tranf.f	5.1  11/6/94
      subroutine tranf(nparm,nz,natoms,ianz,fx,f,ib,b,g,ll)
c***begin prologue     tranf.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)tranf.f	5.1   11/6/94
c***purpose            
c***description
c     routine to transform cartesian first derivatives to
c     derivatives over internal coordinates.
c
c     arguments
c
c     nparm  ... number of z-matrix degrees of freedom (3*nz-6).
c     nz     ... number of rows in the z-matrix.
c     ianz   ... integer atomic numbers of z-matrix elements.
c     fx     ... input vector (length 3*natoms) containing
c                cartesian derivatives.
c     f      ... output vector of length nparm containing
c                derivatives over internal coordinates.
c     ib     ... integer b-matrix as produced by formbg.
c     b      ... b-matrix as produced by formbg.
c     g      ... g-matrix produced by formbg.
c     ll     ... scratch vector of length nz.
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       tranf.f
      implicit none
c     --- input variables -----
      integer nparm,nz,natoms
c     --- input arrays (unmodified) ---
      integer ianz(nz)
      real*8 fx(3*natoms)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      integer ib(4,nparm)
      real*8 f(nparm),b(3,4,nparm),g(nparm,nparm)
c     --- output variables ---
c     --- scratch arrays ---
      integer ll(nz)
c     --- local variables ---
      integer i,j,k,l,ii,k1
      real*8 zero,r
      parameter (zero=0.0d+00)
c
c
      call rzero(f,nparm)
      ii=1
      do 10 i=1,nz
         ll(i)=0
         if(ianz(i).ne.-1) then
            ll(i)=ii
            ii=ii+1
         endif
   10 continue
c
c
      do 60 i=1,nparm
         r=zero
         do 30 k1=1,4
            k=ib(k1,i)
            if(k.eq.0) goto 40
            k=ll(k)
            if(k.ne.0) then
               k=3*(k-1)
               do 20 l=1,3
                  r=r+b(l,k1,i)*fx(k+l)
   20          continue
            endif
   30    continue
c
   40    do 50 j=1,nparm
            f(j)=f(j)+g(i,j)*r
   50    continue
   60 continue
c
c
      return
      end
