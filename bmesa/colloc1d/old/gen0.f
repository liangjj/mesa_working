c $Header: gen0.f,v 1.2 92/12/12 09:34:07 bis Exp $
*deck gen0
      subroutine gen0 (diag,sudiag,spdiag,rhs,guess,m,type)
      implicit integer(a-z)
      real*8 diag, sudiag, spdiag, rhs, guess, invd
      character*(*) type
      dimension diag(0:m), sudiag(0:m), spdiag(0:m)
      dimension rhs(0:m), guess(0:m)
      common /io/ inp, iout
      write(iout,*) '     guess type = ',type
*
*          subroutine to generate a starting guess for the matrix iteration
*          based on pre-conditioning.
*
*          this subroutine is specifically designed for a 
*          tridiagonal system.
*
      rhs(0)=0.d0 
      if (type.eq.'unit-matrix') then
          call copy(rhs(1),guess(1),m)
      elseif (type.eq.'diagonal') then
          do 10 i=1,m
             invd=1.d0/diag(i)
             rhs(i)=rhs(i)*invd
   10     continue
          call copy(rhs(1),guess(1),m)
      elseif (type.eq.'lower-triangular') then
          rhs(1)=rhs(1)/diag(1) 
          do 20 i=2,m
             rhs(i)=(rhs(i)-sudiag(1)*rhs(i-1))/diag(i)
   20     continue
          call copy(rhs(1),guess(1),m)
      else
          call lnkerr('you screwed up in guess routine')
      endif
      return
      end


