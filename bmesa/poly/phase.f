*deck phase.f
c***begin prologue     phase
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           phase shift for s-waves
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***references         
c
c***routines called    
c***end prologue       phase
      subroutine phase(srf,eig,energy,rbox,nrt)
      implicit integer (a-z)
      real*8 srf, eig, energy, rbox, k, pshft, rmat, pi
      dimension srf(nrt), eig(nrt)
      common/io/inp, iout
      data pi/3.141592653d0/
      rmat=0.d0
      do 10 i=1,nrt
         rmat = rmat + srf(i)*srf(i)/( eig(i) -energy )
   10 continue
      rmat=.5d0*rmat
      write(iout,1) energy, rmat
      k=sqrt(2.d0*energy)
      pshft=atan(k*rmat) - k*rbox
      if(pshft.lt.0.d0) then
         pshft=pi+pshft
      endif          
      write(iout,2) energy, pshft
      return
    1 format(/,5x,'energy = ',e15.8,2x,'r-matrix = ',e15.8)      
    2 format(/,5x,'energy = ',e15.8,2x,'phase shift = ',e15.8)      
      end       
