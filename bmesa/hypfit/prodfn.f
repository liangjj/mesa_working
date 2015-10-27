*deck prodfn.f
c***begin prologue     prodfn
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           product of hyperspherical angular and radial
c***                   function and first derivative.
c***author             schneider, barry (nsf)
c***source             
c***purpose            compute product of angular and radial
c***                   hyperspherical functions and its first
c***                   derivative wrt hyperspherical variables
c***                   from the individual functions and their derivatives.
c***description        
c***                   fp =  output product functions
c***                   dfp1 = output derivative wrt angle
c***                   dfp2 = output derivative wrt radius
c***                   fang = input angular function  
c***                   dfang = derivative of angular function  
c***                   frad = input radial function  
c***                   dfrad = input derivative of radial function  
c***references       
c
c***routines called
c***end prologue       prodfn
      subroutine prodfn(fp,dfp1,dfp2,fang,dfang,frad,dfrad,n1,n2)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 fp, dfp1, dfp2, fang, dfang, frad, dfrad
      dimension fp(n2,n1), dfp1(n2,n1), dfp2(n2,n1)
      dimension fang(n2,n1), dfang(n2,n1), frad(n2,n1), dfrad(n2,n1)
      call vmul(fp,fang,frad,n1*n2)
      call vmul(dfp1,frad,dfang,n1*n2)
      call vmul(dfp2,dfrad,fang,n1*n2)
      return
      end
