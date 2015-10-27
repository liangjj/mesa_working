*deck init.f
c***begin prologue     init
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           non-linear schroedinger equation, gross-pitaevski
c***author             schneider, barry (nsf)
c***source             timeprp
c***purpose            starting approximation to gp equation
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       init
      subroutine init(ham0,ham,eig,tmat,psi0,scale,n)
      implicit integer (a-z)
      real*8 ham0, ham, eig, psi0, tmat, scale
      dimension ham0(n,n), ham(n,n), eig(n), psi0(n), tmat(n,n)
      common/io/inp, iout
      call copy(ham0,ham,n*n)
      call tred2(n,n,ham,eig,psi0,tmat)
      call tql2(n,n,eig,psi0,tmat,ierr)
      call copy(tmat(1,1),psi0,n)
      write(iout,1) eig(1)/scale
      return
 1    format(/,5x,'initializing the wavefunction. energy = ',e15.8)
      end       
