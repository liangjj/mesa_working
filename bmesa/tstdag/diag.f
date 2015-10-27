*deck diag.f
c***begin prologue     diag
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           diagonalize
c***author             schneider, barry (nsf)
c***source             math
c***purpose            diagonaliztion driver
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       diag
      subroutine diag(ham,eig,work,dum,n,typ)
      implicit integer (a-z)
      real*8 ham, eig, work, dum
      character*80 title
      dimension ham(n,n), eig(n), work(*), dum(*)
      common/io/inp, iout 
      call tred2(n,n,ham,eig,work,ham)
      call tql2(n,n,eig,work,ham,ierr)
      title='eigenvalues'
      call prntrm(title,eig,n,1,n,1,iout)
      return
      end       

