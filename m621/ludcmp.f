*deck @(#)ludcmp.f	5.1  11/6/94
C  (C) Copr. 1986-92 Numerical Recipes Software 4.
      SUBROUTINE ludcmp(a,n,np,indx,d)
      implicit none

c....	given a matrix a(1:n,1:n), with physical dimension np by np,
c	this routine replaces it by the LU decomposition of a row-
c	wise permutation of itself.  a and n are input.  a is output,
c	arranged a s in equation (2.3.14) above;  indx(1:n) is an output
c	vector that records the row permutation effected by the partial
c	pivoting;  d is output as +/- 1 depending on whether the
c	number of row interchanges was even or odd, respectively.  this
c	routine is used in combination with DLUBKSB to solve linear
c	equations or invert a matrix
      INTEGER n,np,indx(n),NMAX
      REAL*8 d,a(np,np),TINY
      real*8 zero,one
      parameter (zero=0.0d+00,one=1.0d+00)
      PARAMETER (NMAX=4096,TINY=1.0d-20)

C	LARGEST EXPECTED N, AND A SMALL NUMBER

      INTEGER i,imax,j,k
      REAL*8 aamax,dum,sum,vv(NMAX)

C	VV STORES THE IMPLICIT SCALING OF EACH ROW

      d=one

C	NO ROW INTERCHANGES YET

      do 5 i=1,n
        aamax=zero

C	LOOP OVER ROWS TO GET IMPLICIT SCALING INFORMATION

        do 10 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
  10    continue

        if (aamax.eq.zero) then
           call lnkerr('singular matrix in ludcmp')
        endif
        vv(i)=one/aamax
  5   continue

C	THIS IS THE LOOP OVER COLUMNS OF CROUT'S METHOD

      do 15 j=1,n

C	THIS IS EQUATION (2.3.12) EXCEPT FOR I=J

        do 20 i=1,j-1
          sum=a(i,j)

          do 25 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
  25      continue

          a(i,j)=sum
  20    continue

C	INITIALIZE FOR THE SEARCH FOR LARGEST PIVOT ELEMENT

        aamax=zero

C	THIS IS I=J OF EQUATION (2.3.12) AND I-J+1,...,N OF EQUATION
C	(2.3.13)

        do 30 i=j,n
          sum=a(i,j)

          do 35 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
  35      continue

C	FIGURE OF MERIT FOR THE PIVOT

          a(i,j)=sum
          dum=vv(i)*abs(sum)

C	IS IT BETTER THAN THE BEST SO FAR?

          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
  30    continue

C	DO WE NEED TO INTERCHANGE ROWS?

        if (j.ne.imax)then

C	YES, DO SO

          do  40 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
  40      continue

C	AND CHANGE THE PARITY OF D

          d=-d

C	ALSO CHANGE THE SCALE FACTOR

          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.zero)a(j,j)=TINY

C	IF THE PIVOT ELEMENT IS ZERO THE MARTIX IS SINGULAR (AR LEAST TO
C	THE PRECISION OF THE ALGORITHM).  FOR SOME APPLICATIONS OF
C	SINGULAR MATRICES, IT IS DESIRABLE TO SUBSTITUTE TINY FOR ZERO.

C	NOW, FINALLY, DIVIDE BY THE PIVOT ELEMENT

        if(j.ne.n)then
          dum=one/a(j,j)

          do 45 i=j+1,n
            a(i,j)=a(i,j)*dum
  45      continue

        endif

  15   continue

      return
      END
