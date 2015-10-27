*deck mktrid.f
c***begin prologue     mktrid
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           lanczos
c***author             schneider, barry (nsf)
c***source             math
c***purpose            lanczos diagonalization
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       mktrid
      subroutine mktrid(ham,eig,work,dum,temp,n,niter,prnt)
      implicit integer (a-z)
      real*8 ham, eig, work, dum, temp
      character*80 title
      logical prnt
      dimension ham(n,n), eig(n), work(*), dum(*), temp(*)
      common/io/inp, iout 
      call lancz(temp,ham,eig,dum,work,n,niter,.true.)
      title='eigenvalues'
      call prntrm(title,eig(2),niter,1,niter,1,iout)
      return
      end       

