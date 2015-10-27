*deck @(#)minv.f	1.1  11/30/90
      subroutine minv(ab,n,nd,scratch,det,eps,m,mode)
c
      implicit integer (a-z)
c
      real*8 ab(nd,n)
      real*8 scratch(n,2)
      real*8 det
      real*8 eps
      real*8 rcond
c
      if (mode.eq.1) then
         call sgeco(ab,nd,n,scratch(1,1),rcond,scratch(1,2))
         if (rcond.eq.0.0d+00) then
            call lnkerr('in minv, the matrix is singular')
         end if
         call sgedi(ab,nd,n,scratch(1,1),det,scratch(1,2),11)
      else
         call lnkerr('unimplemented call to minv')
      end if
c
c
      return
      end
