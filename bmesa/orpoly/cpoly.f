*deck cpoly.f
c***begin prologue     cpoly
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           orthogonal polynomials
c***author             schneider, barry (nsf)
c***source             math
c***purpose            compute coefficients of three term recursion
c***                   relation, points and weights of orthogonal
c***                   polynomials on (-1,1), (0,infinity),
c***                   or (-infinity,infinity).
c***description        the recursion coefficients are computed, and then
c***                   the roots and weights determined by diagonalizing
c***                   the appropriate tridiagonal matrix.
c***                   the recursion relation used is,
c***
c***                   b p (x) = (x - a ) p   (x) - b   p   (x)
c***                    j j            j   j-1       j-1 j-2
c***                   
c***                   where the p's are normalized and orthogonal.
c***references         
c
c***routines called    class, gbtql2, gbslve
c***end prologue       cpoly
      subroutine cpoly(alpha,beta,endpts,rbeg,rend,kpts,x,wts,a,b,
     1                 norm0,scr,n,type,prnt)
      implicit integer (a-z)
      real*8 alpha, beta, endpts, rbeg, rend, x, wts, a, b, scr, pi
      real*8 h, norm0, gbslve, gam, t1
      character*(*) type
      character*80 title
      logical prnt
      dimension x(n), wts(n), a(n), b(n), scr(n)
      dimension endpts(2)
      common/io/inp, iout 
      data pi/3.141592653589793238462643d0/
      if (type.eq.'simpson') then
          kind=0
      elseif(type.eq.'legendre') then
          kind=1
      elseif(type.eq.'chebyshev-1') then
          kind=2
      elseif(type.eq.'chebyshev-2') then
          kind=3
      elseif(type.eq.'hermite') then
          kind=4
      elseif(type.eq.'jacobi') then
          kind=5
      elseif(type.eq.'laguerre') then
          kind=6
      else
          call lnkerr('error in type polynomial')
      endif
      if(kind.eq.0) then
c           this is for simpson rule points and weights
         if(2*(n/2).eq.n) then
            call lnkerr('n must be odd for simpson rule')
         endif
         if(n.le.1) then
            x(1) = 0.d+00
            wts(1) = 2.d+00
         else
            h = 2.d+00/(n-1)
            x(1) = -1.d+00
            x(n) = 1.d+00
            wts(1) = h/3.d+00
            wts(n) = h/3.d+00
            nm1 = n-1
            do 10 i=2,nm1
               x(i) = x(i-1) + h
               wts(i) = 4.d+00 - 2.d+00*(i-2*(i/2))
               wts(i) = wts(i)*h/3.d+00
 10         continue   
         endif
      else
c
c           find the coefficients a and b of the three term recursion
c           relationship.
c
c           the matrix of coefficients is assumed to be symmetric.
c           the array a contains the diagonal elements, the array
c           b the off-diagonal elements.
c           make appropriate changes in the lower right 2 by 2
c           submatrix.
c
         call class (kind, n, alpha, beta, b, a, norm0)
c
         if (kpts.ne.0)  then
             if (kpts.eq.1) then
c
c                kpts=1, only a(n) must be changed
c
                 a(n) =gbslve(endpts(1), n, a, b)*b(n-1)**2 + endpts(1)
             endif
             if (kpts.eq.2) then
c
c                kpts=2, a(n) and b(n-1) must be recomputed
c
                 gam =gbslve(endpts(1), n, a, b)
                 t1 = ((endpts(1) - endpts(2))/
     1                        (gbslve(endpts(2), n, a, b) - gam))
                 b(n-1) =  sqrt(t1)
                 a(n) = endpts(1) + gam*t1
             endif
         endif
c
c           note that the indices of the elements of b run from 1 to n-1
c           and thus the value of b(n) is arbitrary.
c           now compute the eigenvalues of the symmetric tridiagonal
c           matrix, which has been modified as necessary.
c           the method used is a ql-type method with origin shifting
c
         wts(1) = 1.0d0
         do 40 i = 2, n
            wts(i) = 0.0d0
 40      continue   
c
         call copy(a,x,n)
         call copy(b,scr,n)
         call gbtql2 (n, x, scr, wts, ierr)
         do 50 i = 1, n
            wts(i) = norm0 * wts(i) * wts(i)
 50      continue   
      endif
c
      call iosys('write integer "number of polynomial points" to '//
     1           'lamdat',1,n,0,' ')
      call iosys('write real "polynomial points" to lamdat',n,
     1            x,0,' ')
      call iosys('write real "polynomial weights" to lamdat',n,
     1           wts,0,' ')
      if (kind.ne.0) then
          call iosys('write real "polynomial a coefficients" to lamdat',
     1                n,a,0,' ')
          call iosys('write real "polynomial b coefficients" to lamdat',
     1                n,b,0,' ')
      endif
c
      norm0=1.d0/sqrt(norm0)
      if(prnt) then
          write(iout,1) norm0
          title='polynomial points'
          call prntrm(title,x,n,1,n,1,iout)
          title='polynomial weights'
          call prntrm(title,wts,n,1,n,1,iout)
          if (kind.ne.0) then
              title='a coefficients'
              call prntrm(title,a,n,1,n,1,iout)
              title='b coefficients'
              call prntrm(title,b,n-1,1,n-1,1,iout)
          endif
      endif
      return
 1    format(/,'normalization constant for p(0) = ',e15.8)
      end       
