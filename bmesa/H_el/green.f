*deck green.f 
c***begin prologue     green
c***date written       020303   (yymmdd)
c***revision date               (yymmdd)
c***keywords           
c***                   
c***author             schneider, b. i.(nsf)
c***source             hcoll
c***purpose            greens function elements for s-waves
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       green
      subroutine green(r_k,i_k,r,ene,n)
c
      implicit integer (a-z)
      real*8 r_k, i_k, ene, r
      real*8 k, fac
      dimension r(n), r_k(n), i_k(n)
      common/io/inp, iout      
      k=sqrt(2.d0*ene)
      fac = sqrt(2.d0/k)
      do 10 i=1,n
         r_k(i) = fac*sin(k*r(i))
         i_k(i) = fac*cos(k*r(i))
 10   continue   
      return
      end


















