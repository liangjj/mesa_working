*deck itham.f 
c***begin prologue     itham
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           
c***                   hamiltonian, dvr, hyperspherical
c***author             schneider, b. i.(nsf)
c***source             
c***purpose            explicit construction of total hyperspherical
c***                   hamiltonian in dvr basis.
c***                   
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       itham
      subroutine itham(ham,vad,rho,n)
c
      implicit integer (a-z)
      real*8 ham, vad, rho, two, rhosq 
      dimension ham(n,n), vad(n)
      common/io/inp, iout      
      data two / 2.d0 /
      rhosq=1.d0/(rho*rho)
      do 10 i=1,n
         ham(i,i) = ( ham(i,i) - two )*rhosq + vad(i)    
 10   continue   
c
      return
      end
