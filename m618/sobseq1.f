*deck @(#)sobseq1.f	5.1  11/6/94
      SUBROUTINE sobseq1(n,x)
      implicit none
      INTEGER n,MAXBIT,MAXDIM
      REAL*8 x(*)
      PARAMETER (MAXBIT=30,MAXDIM=12)

C...	WHEN N IS NEGATIVE, INTERNALLY INITIALIZES A SET OF MAXBIT
C	DIRECTION NUMBERS FOR EACH OF MAXDIM DIFFERENT SOBOL' SEQUENCES.
C	WHEN N IS POSITIVE (BUT <= MAXDIM), RETURNS AS THE VECTOR
C	X(1..N) THE NEXT VALUES FROM N OF THESE SEQUENCES.  (N MUST
C	NOT BE CHANGES BETWEEN INITIALIZATIONS.)

      INTEGER i,im,in,ipp,j,k,l,ip(MAXDIM),iu(MAXDIM,MAXBIT)
     &,	iv(MAXBIT*MAXDIM),ix(MAXDIM),mdeg(MAXDIM),xor
      REAL*8 fac
      real*8 one,two
      parameter (one=1.0d+00,two=2.0d+00)
      SAVE ip,mdeg,ix,iv,in,fac

C...	TO ALLOW BOTH 1D AND 2D ADDRESSING

      EQUIVALENCE (iv,iu)
      DATA 	ip /0,1,1,2,1,4,2,4,7,11,13,14/
     &, 	mdeg /1,2,3,3,4,4,5,5,5,5,5,5/

      if (n.lt.0) then

C...	ALTERNATIVE INITIALIZATION SO THAT REPEATED CALLS WITH
C	THE SAME ARGUMENTS GIVE THE SAME ANSWER

	do 5 i = 1, MAXDIM
		ix(i) = 0
  5     continue

	do 7 i = 1, MAXBIT*MAXDIM
		iv(i)=0
  7     continue

	iu(1,1) = 1
	iu(1,2) = 3
	iu(1,3) = 5
	iu(1,4) = 15
	iu(1,5) = 17

	iu(2,1) = 1
	iu(2,2) = 1
	iu(2,3) = 7
	iu(2,4) = 11
	iu(2,5) = 13

	iu(3,1) = 1
	iu(3,2) = 3
	iu(3,3) = 7
	iu(3,4) = 5
	iu(3,5) = 7

	iu(4,1) = 1
	iu(4,2) = 3
	iu(4,3) = 3
	iu(4,4) = 15
	iu(4,5) = 5

	iu(5,1) = 1
	iu(5,2) = 1
	iu(5,3) = 3
	iu(5,4) = 13
	iu(5,5) = 25

	iu(6,1) = 1
	iu(6,2) = 1
	iu(6,3) = 5
	iu(6,4) = 9
	iu(6,5) = 3

	iu(7,1) = 1
	iu(7,2) = 1
	iu(7,3) = 1
	iu(7,4) = 1
	iu(7,5) = 13

	iu(8,1) = 1
	iu(8,2) = 1
	iu(8,3) = 1
	iu(8,4) = 3
	iu(8,5) = 17

	iu(9,1) = 1
	iu(9,2) = 1
	iu(9,3) = 1
	iu(9,4) = 5
	iu(9,5) = 19

	iu(10,1) = 1
	iu(10,2) = 1
	iu(10,3) = 3
	iu(10,4) = 7
	iu(10,5) = 23

	iu(11,1) = 1
	iu(11,2) = 1
	iu(11,3) = 5
	iu(11,4) = 11
	iu(11,5) = 29

	iu(12,1) = 1
	iu(12,2) = 3
	iu(12,3) = 7
	iu(12,4) = 13
	iu(12,5) = 31

C...	STORED VALUES ONLY REQUIRE NORMALIZATION

        do 9 k=1,MAXDIM

          do 11 j=1,mdeg(k)
            iu(k,j)=iu(k,j)*2**(MAXBIT-j)
  11      continue

C...	USE THE RECURRENCE TO GET OTHER VALUES.

          do 13 j=mdeg(k)+1,MAXBIT
            ipp=ip(k)
            i=iu(k,j-mdeg(k))
            i=xor(i,i/2**mdeg(k))

            do 14 l=mdeg(k)-1,1,-1
              if(iand(ipp,1).ne.0)i=xor(i,iu(k,j-l))
              ipp=ipp/2
  14        continue
            iu(k,j)=i
  13      continue
  9     continue
        fac=one / two**MAXBIT
        in=0

      else
C...	CALCULATE THE NEXT VECTOR IN THE SEQUENCE

        im=in

        do 15 j=1,MAXBIT

C...	FIND THE RIGHTMOST ZERO BIT

          if(iand(im,1).eq.0)goto 1
          im=im/2
15      continue
        call lnkerr('MAXBIT too small in sobseq1')
1       im=(j-1)*MAXDIM

C...	XOR THE APPROPRIATE DIRECTION NUMBER INTO EACH COMPONENT OF THE
C	VECTOR AND CONVERT TO A FLOATING NUMBER

        do 16 k=1,min(n,MAXDIM)
          ix(k)=xor(ix(k),iv(im+k))
          x(k)=ix(k)*fac
16      continue

C...	INCREMENT THE COUNTER

        in=in+1

      endif
c
c
      return
      END
