*deck potnl.f
c***begin prologue     potnl
c***date written       970813   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           non-linear schroedinger equation, gross-pitaevski
c***author             schneider, barry (nsf)
c***source             trap3d
c***purpose            form current approximation to non-linear
c***                   potential.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       potnl
      subroutine potnl(vnl,eig,pndvr,psiin,gamma,n)
      implicit integer (a-z)
      real*8 vnl, eig, pndvr, psiin
      real*8 gamma
      dimension vnl(n), eig(n), pndvr(n,0:n-1)
      dimension psiin(n)
      common/io/inp, iout
c      write(iout,*) vnl
c      write(iout,*) eig
c      write(iout,*) pndvr
c      write(iout,*) gamma

c     transform wavefunction
c
      call ebc(vnl,pndvr,psiin,n,n,1)
c
c     calculate non-linear term
c
      do 10 i=1,n
         vnl(i) = vnl(i)/eig(i)
 10   continue
      do 20 i=1,n
         vnl(i) = vnl(i)*vnl(i)
 20   continue   
      call smul(vnl,vnl,gamma,n) 
      return 
      end       
