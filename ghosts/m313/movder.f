*deck @(#)movder.f	1.1  11/30/90
      subroutine movder(der,grad,nat,dercen,npass)
c
      implicit integer (a-z)
c
      real*8 der(3,4),grad(3,nat)
      integer dercen(4)
c
      if (npass.eq.0) return
c
      do 10 c=1,3
         if (npass.eq.1) then
            grad(c,dercen(1))=grad(c,dercen(1))+der(c,1)
            grad(c,dercen(2))=grad(c,dercen(2))+der(c,2)
            grad(c,dercen(3))=grad(c,dercen(3))+der(c,3)
            grad(c,dercen(4))=grad(c,dercen(4))-der(c,1)-der(c,2)-
     #                                          der(c,3)
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
