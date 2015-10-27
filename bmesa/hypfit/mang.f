*deck mang.f
c***begin prologue     mang
c***date written       001230   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            hyperspherical angular functions.
c***                   
c***                   Y(alpha,m) = N*sin[ (2*m+2)*alpha ]
c***                                      / { sin[alpha]*cos[alpha] }
c***                      N = 2/sqrt(pi) , m = integer
c***
c***                   Y(alpha,m) = 2*N*sin[ (2*m+2)*alpha ]
c***                                      / sin[2*alpha]
c***references         
c
c***routines called    
c***end prologue       mang
      subroutine mang(yalf,dyalf,alf,nrmalf,pi2,fac,m)
      implicit integer (a-z)
      real*8 yalf, dyalf, alf, nrmalf
      real*8 pi2, fac, angle, sialf, calf, smalf, cmalf
      logical prn
      common/io/inp, iout
      mbeg=m+2
      imul=m-2*(m/2)
      if(imul.eq.0) then
         imul=1
      else
         imul=-1
      endif 
      if(alf.eq.0.d0) then
c
c        by l'hospital rule 
c
         yalf = nrmalf*mbeg
         dyalf = 0.d0
      elseif(alf.eq.pi2) then   
c
c        by l'hospital rule 
c
         yalf = imul*nrmalf*mbeg
         dyalf = 0.d0
      else   
         angle=mbeg*alf
         sialf=sin(2.d0*alf)
         calf=cos(2.d0*alf)
         smalf=sin(angle)
         cmalf=cos(angle)
         yalf = fac*smalf/sialf
         dyalf = fac*(mbeg*cmalf*sialf 
     1                    - 
     2           2*smalf*calf)/(sialf*sialf)
      endif
      return
      end       
