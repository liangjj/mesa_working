*deck @(#)rveceq.f	5.1  11/6/94
      function rveceq(a,b,n)
c***begin prologue     rveceq.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             unknown
c***source             @(#)rveceq.f	5.1   11/6/94
c***purpose            compares two real*8 vectors 
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       rveceq.f
      implicit none
c     --- input variables -----
      integer n
c     --- input arrays (unmodified) ---
      real*8 a(n),b(n)
c     --- input arrays (scratch) ---
c     --- function returned ---
      logical rveceq
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i
c
      rveceq=.true.
      do 1 i=1,n
         if (a(i).ne.b(i)) then
            rveceq=.false.
            return
         endif
    1 continue
c
c
      return
      end
