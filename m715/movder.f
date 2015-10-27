*deck @(#)movder.f	5.1  11/28/95
      subroutine movder(der,grad,nat,dercen,npass)
c***begin prologue     movder.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             saxe, paul(lanl)
c***source             @(#)movder.f	5.1   11/28/95
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       movder.f
      implicit none
c     --- input variables -----
      integer nat,npass
c     --- input arrays (unmodified) ---
      integer dercen(4)
      real*8 der(3,4)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 grad(3,nat)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer c
c
      if (npass.eq.0) return
c
      do 10 c=1,3
         if (npass.eq.1) then
            grad(c,dercen(1))=grad(c,dercen(1))+der(c,1)
            grad(c,dercen(2))=grad(c,dercen(2))+der(c,2)
            grad(c,dercen(3))=grad(c,dercen(3))+der(c,3)
            grad(c,dercen(4))=grad(c,dercen(4))-der(c,1)-der(c,2)-
     $                                          der(c,3)
         else if (npass.eq.2) then
            grad(c,dercen(1))=grad(c,dercen(1))+der(c,1)
            grad(c,dercen(2))=grad(c,dercen(2))+der(c,2)
            grad(c,dercen(3))=grad(c,dercen(3))-der(c,1)-der(c,2)
         else if (npass.eq.3) then
            grad(c,dercen(1))=grad(c,dercen(1))+der(c,1)
            grad(c,dercen(3))=grad(c,dercen(3))+der(c,2)
            grad(c,dercen(2))=grad(c,dercen(2))-der(c,1)-der(c,2)
         else if (npass.eq.4) then
            grad(c,dercen(1))=grad(c,dercen(1))+der(c,1)
            grad(c,dercen(2))=grad(c,dercen(2))-der(c,1)
         end if
   10 continue
c
c
      return
      end
