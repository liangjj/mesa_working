*deck diag
      subroutine diag(ham,eig,scr,ncon)
      implicit integer (a-z)
      real*8 ham, eig, scr
      dimension ham(ncon,ncon), eig(ncon), scr(*)
      common /io/ inp, iout
      call tred2(ncon,ncon,ham,eig,scr,ham)
      call tql2(ncon,ncon,eig,scr,ham,ierr)
      write(iout,1) (eig(i),i=1,ncon)  
 1    format(/,1x,'r matrix eigenvalues',(/,5e15.8))
      return
      end
