*deck guess.f
      subroutine guess(ham,eig,work,n,m,nroots)
      implicit integer (a-z)
      real*8 ham, eig, work
      character*80 title
      dimension ham(n,*), eig(*), work(*)
      common/io/inp, iout 
      call tred2(n,m,ham,eig,work,ham)
      call tql2(n,m,eig,work,ham,ierr)
      title='guess eigenvalues'
      call prntrm(title,eig,nroots,1,nroots,1,iout)
      return
      end       

