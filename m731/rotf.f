*deck @(#)rotf.f	5.1  11/6/94
      subroutine rotf(natoms,tr,fin,fout)
c***begin prologue     rotf.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkely, et al. (g82)
c***source             @(#)rotf.f	5.1   11/6/94
c***purpose            
c***description
c     routine to rotate energy derivatives back to original
c     cartesian axes.  note that this is a general routine,
c     and is not restricted to derivatives.
c
c     arguments:
c
c     natoms ... number of atoms.
c     tr     ... (3 by 3) rotation matrix.
c     fin    ... input vector of length 3*natoms,
c     fout   ... output vector, can be same as input vector.
c     
c
c***references
c
c***routines called
c
c***end prologue       rotf.f
      implicit none
c     --- input variables -----
      integer natoms
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
      real*8 fin(3,natoms)
c     --- output arrays ---
      real*8 fout(3,natoms)
c     --- output variables ---
c     --- scratch arrays ---
      real*8 tr(3,3)
c     --- local variables ---
      integer iat,i
      real*8 tx(3)
c
c
      do 30 iat=1,natoms
         do 10 i=1,3
            tx(i)=tr(1,i)*fin(1,iat) + tr(2,i)*fin(2,iat) +
     $            tr(3,i)*fin(3,iat)
   10    continue
         do 20 i=1,3
            fout(i,iat)=tx(i)
   20    continue
   30 continue
c
c
      return
      end
