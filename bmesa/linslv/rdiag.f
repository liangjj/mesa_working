*deck rdiag.f
c***begin prologue     rdiag
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           diagonalization
c***author             schneider, barry (nsf)
c***source             
c***purpose            driver for real diagonalization.
c***                   
c***references         
c
c***routines called    
c***end prologue       rdiag
      subroutine rdiag(ham,u,eig,work,n)
      implicit integer (a-z)
      real*8 ham, u, eig, work
      character*80 title
      dimension ham(n,n), u(n,n), eig(n), work(*)
      common/io/inp, iout 
      call copy(ham,u,n*n)
      call tred2(n,n,u,eig,work,u)
      call tql2(n,n,eig,work,u,ierr)
      title='eigenvalues of zeroth order hamiltonian'
      call prntrm(title,eig,n,1,n,1,iout)
      return
      end       
