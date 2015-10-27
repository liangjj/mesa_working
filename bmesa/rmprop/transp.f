      subroutine transp (a,b,n,nmax)
c
c*    subroutine transp( a, b, n, nmax )
c**   _transp" -- copies the transpose of a(nmax,n) to b(nmax,n)
c**   a(nmax,n) -- matrix to be transposed  (real)
c**   b(nmax,n) -- output matrix, the transpose of a(nmax,n)  (real)
c**   nmax, n -- apropriate dimension variables   (integer*4)
c
      real*8 a, b
      dimension a(nmax,n), b(nmax,n)
      do 10 i=1,n
      do 10 j=1,n
      b(j,i)=a(i,j)
   10 continue
      return
      end
