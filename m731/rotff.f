*deck @(#)rotff.f	5.1  11/6/94
      subroutine rotff(natoms,tr,ffxin,ffxout)
c***begin prologue     rotff.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)rotff.f	5.1   11/6/94
c***purpose            
c***description
c     routine to rotate the second derivative matrix (stored
c     in lower triangular form) from one axis system to another.
c
c     arguments:
c
c     natoms ... number of atoms.
c     tr     ... (3 by 3) rotation matrix.
c     ffxin  ... input array of second-derivatives.
c     ffxout ... output array, can be same as input array.
c
c***references
c
c***routines called
c
c***end prologue       rotff.f
      implicit none
c     --- input variables -----
      integer natoms
c     --- input arrays (unmodified) ---
      real*8 ffxin(natoms*(natoms+1)/2)
c     --- input arrays (scratch) ---
      real*8 tr(3,3)
c     --- output arrays ---
      real*8 ffxout(natoms*(natoms+1)/2)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer lind
      integer iat,jat,i,j,k,l,m,locat
      real*8 t(3,3), t1(3,3)
      real*8 zero
c
      parameter (zero=0.0d+00)
c
      lind(i,j)=((max(i,j)*(max(i,j)-1)/2)+min(i,j))
c
c     --- loop over all atoms.
      do 50 iat=1,natoms
         i=3*(iat-1)
         do 50 jat=1,iat
            j=3*(jat-1)
c
c           pluck out the current (3*3) matrix.
            do 10 k=1,3
               do 10 l=1,3
                  locat=lind(k+i,l+j)
   10             t(k,l)=ffxin(locat)
c
c           transform by tr.
            do 20 k=1,3
               do 20 l=1,3
                  t1(k,l)=zero
                     do 20 m=1,3
   20                   t1(k,l)=t1(k,l) +tr(m,l)*t(k,m)
            do 30 l=1,3
               do 30 k=1,3
                  t(k,l)=zero
                     do 30 m=1,3
   30                   t(k,l)=t(k,l) +tr(m,k)*t1(m,l)
c
c           pack the matrix back into ffx.
            do 40 k=1,3
               do 40 l=1,3
                  locat=lind(k+i,l+j)
   40             ffxout(locat)=t(k,l)
   50    continue
c
c
      return
      end
