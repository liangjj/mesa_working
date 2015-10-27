*deck @(#)lubksb.f	5.1  11/6/94
C  (C) Copr. 1986-92 Numerical Recipes Software 4.
      SUBROUTINE lubksb(a,n,np,indx,b)
      implicit none

C...	SOLVES THE SET OF N LINEAR EQUATIONS AX=B.  HERE A IS INPUT,
C	NOT AS THE MATRIX A BUT RATHER AS ITS LU DECOMPOSITION,
C	DETERMINED BY THE ROUTINE DLUDCMP.  INDX IS INPUT AS THE
C	PERMUTATION VECTOR RETURNED BY DLUDCMP.  B(1:N) IS INPUT AS THE
C	RIGHT SIDE VECTOR B, AND RETURNS WITH THE SOLUTION VECTOR X.
C	A, N, NP, AND INDX ARE NOT MODIFIED BY THIS ROUTINE AND CAN BE
C	LEFT IN PLACE FOR SUCCESSIVE CALLS WITH DIFFERENT RIGHT SIDES B.
C	THIS ROUTINE TAKES INTO ACCOUNT THE POSSIBILITY THAT B WILL
C	BEGIN WITH MANY ZERO ELEMENTS, SO IT IS EFFICIENT FOR USE IN
C	MATRIX INVERSION.
      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL*8 sum
      real*8 zero
      parameter (zero=0.0d+00)

C	WHEN II IS SET TO A POSITIVE VALUE, IT WILL BECOME THE INDEX OF
C	THE FIRST NONVANISHING ELEMENT OF B.  WE NOW DO THE FORWARD
C	SUBSTITUTION, EQUATION (2.3.6).  THE ONLY NEW WRINKLE IS TO
C	UNSCRAMBLE THE PERMUTATION AS WE GO.

      ii=0
      do 5 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then

          do 10 j=ii,i-1
            sum=sum-a(i,j)*b(j)
  10      continue

        else if (sum.ne.zero) then

C	A NONZERO ELEMENT WAS ENCOUNTERED, SO FROM NOW ON WE WILL HAVE
C	TO DO THE SUMS IN THE LOOP ABOVE.

          ii=i
        endif
        b(i)=sum
  5   continue

C	NOW WE DO THE BACKSUBSTITUTION, EQUATION (2.3.7)

      do 15 i=n,1,-1
        sum=b(i)
	
        do 20 j=i+1,n
          sum=sum-a(i,j)*b(j)
  20    continue

C	STORE A COMPONENT OF THE SOLUTION VECTOR X

        b(i)=sum/a(i,i)

  15  continue

      return
      END
