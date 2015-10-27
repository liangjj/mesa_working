*deck @(#)runit.f	5.1  11/6/94
      subroutine runit(a,n)
c***begin prologue     runit.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             unknown
c***source             @(#)runit.f	5.1   11/6/94
c***purpose            fill a with the unit matrix. 
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       runit.f
      implicit none
c     --- input variables -----
      integer n
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 a(n,n)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i
c
      call rzero(a,n**2)
      do 1 i=1,n
         a(i,i)=1.0d+00
    1 continue
c
c
      return
      end
