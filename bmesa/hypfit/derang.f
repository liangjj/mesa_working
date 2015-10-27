*deck derang.f
c***begin prologue     derang
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           chain rule for df/dang from df/dx and df/dy
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***description          
c***references       
c
c***routines called
c***end prologue       derang
      subroutine derang(dfdang,dfdr1,dfdr2,rho,ang,n1,n2)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 dfdang, dfdr1, dfdr2, rho, ang
      character*80 title
      dimension dfdang(n2,n1), dfdr1(n2,n1), dfdr2(n2,n1)
      dimension rho(n2,n1), ang(n2,n1)
      do 10 i=1,n1
         do 20 j=1,n2
            dfdang(j,i) = rho(j,i)*( dfdr1(j,i)*cos(ang(j,i)) 
     1                            - 
     2                    dfdr2(j,i)*sin(ang(j,i)) )
 20      continue   
 10   continue
      return
      end
