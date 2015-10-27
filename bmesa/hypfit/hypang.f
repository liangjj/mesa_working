*deck hypang.f
c***begin prologue     hypang
c***date written       001230   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            hyperspherical angular functions in two dimensions
c***                   
c***                   Y(alpha,m) = N*sin[ (2*m+2)*alpha ]
c***                                      / { sin[alpha]*cos[alpha] }
c***                      N = 2/sqrt(pi) , m = integer
c***
c***references         
c
c***routines called    
c***end prologue       hypang
      subroutine hypang(ang,yangm,dyangm,m,n1,n2)
      implicit integer (a-z)
      real*8 ang, yangm, dyangm 
      real*8 pi, nrmalf, pi2, fac, two
      dimension ang(n2,n1), yangm(n2,n1), dyangm(n2,n1)
      common/io/inp, iout
      data pi  / 3.141592653589793238462643d+00 /
      data pi2 / 1.5707963267948966192313215d+00 /
      data two / 2.d0 /
      nrmalf=two/sqrt(pi)
      fac=two*nrmalf
c
      do 10 i=1,n1
         do 20 j=1,n2
c
c           compute the angular functions
c           be careful when the angle is zero or pi/2
c
            call mang(yangm(j,i),dyangm(j,i),ang(j,i),nrmalf,pi2,fac,m)
 20      continue   
 10   continue   
      return
      end       















