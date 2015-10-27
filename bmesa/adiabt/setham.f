*deck setham.f 
c***begin prologue     setham
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           hyperspherical angular hamiltonian
c***author             schneider, b. i.(nsf)
c***source             
c***purpose            using dvr functions set up the
c***                   hyperspherical hamiltonian.
c***                   
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       setham
      subroutine setham(hamad,ham,vad,rho,n)
c
      implicit integer (a-z)
      real*8 hamad, ham, vad, rho, two, rhosq 
      character*80 title
      dimension hamad(n,n), ham(n,n), vad(n)
      common/io/inp, iout      
      data two / 2.d0 /
      rhosq=rho*rho
      call copy(ham,hamad,n*n)
      do 10 i=1,n
         hamad(i,i) = hamad(i,i) - two    
 10   continue
      do 20 i=1,n
         hamad(i,i) = hamad(i,i) + vad(i)*rhosq
 20   continue   
c
      return
      end
