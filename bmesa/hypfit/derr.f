*deck derr.f
c***begin prologue     derr
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           chain rule for df/dr from df/dx and df/dy
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***description          
c***references       
c
c***routines called
c***end prologue       derr
      subroutine derr(dfdr,dfdr1,dfdr2,ang,n1,n2)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 dfdr, dfdr1, dfdr2, ang
      dimension dfdr(n2,n1), dfdr1(n2,n1), dfdr2(n2,n1), ang(n2,n1)
      do 10 i=1,n1
         do 20 j=1,n2
            dfdr(j,i) = dfdr1(j,i)*sin(ang(j,i)) 
     1                         + 
     2                  dfdr2(j,i)*cos(ang(j,i))
 20      continue   
 10   continue
      return
      end
