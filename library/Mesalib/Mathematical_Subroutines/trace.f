*deck @(#)trace.f	5.1  11/6/94
      function trace(a,n)
c***begin prologue     trace.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             unknown
c***source             @(#)trace.f	5.1   11/6/94
c***purpose            traces a matrix
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       trace.f
      implicit none
c     --- input variables -----
      integer n
c     --- input arrays (unmodified) ---
      real*8 a(n,n)
c     --- input arrays (scratch) ---
c     --- function returned ---
      real*8 trace
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i
      real*8 t
c
c
      t=0.0d+00
      do 1 i=1,n
         t=t+a(i,i)
    1 continue
      trace=t
c
      return
      end
