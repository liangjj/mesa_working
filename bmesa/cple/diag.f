*deck diag
c***begin prologue     diag
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           eigenvalues, kinetic
c***author             schneider, barry (nsf)
c***source             
c***purpose            diagonalize hamiltonian matrix
c***description        
c***references       
c
c***routines called
c***end prologue       diag
      subroutine diag(ham,eig,dum,n,prnt)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 ham, eig, dum
      logical prnt
      dimension ham(n,n), eig(n), dum(n)
      call tred2(n,n,ham,eig,dum,ham)
      call tql2(n,n,eig,dum,ham,ierr)
      if (prnt) then
          write(iout,1) (eig(j),j=1,n)
      endif          
      call iosys ('write integer "size of hamiltonian matrix" to '//
     1            'lamdat',1,n,0,' ')
      call iosys ('write real "hamiltonian eigenvalues" to lamdat',
     1             n,eig,0,' ')
    1 format(/,1x,'hamiltonian eigenvalues',
     1         (/,1x,5e15.8))
      return
      end
