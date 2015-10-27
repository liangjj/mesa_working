!=======================================================================
  SUBROUTINE wavefn 
!=======================================================================
!
!       This routine initializes radial functions.  The file 'wfn.inp"
!  is read.  The first radial function with the same electron label is
!  scaled to the current Z as an initial estimate.  Otherwise, an 
!  estimate is formed as a screened hydrogenic function witn ZZ=Z-s(i)
!  and a Hartree-orbital computed, where previously defined orbitals
!  define the potential and the new orbital is constrained to be orthogonal
!  to existing orbitals.
!
!----------------------------------------------------------------------
!
    IMPLICIT NONE
!
      CHARACTER*3 EL1
      CHARACTER*6 AT,TT,ATM(NWD),TRM(NWD)
      CHARACTER*24 TITLE
      LOGICAL LD
!
      nint=ns-ks+1
      call facsb(nt,kx,ks,ns,sb,bs)
!
!  ***** READ THE WAVEFUNCTIONS
!
      IF (IUF .EQ. 0) GO TO 5
2     READ(IUF,END=5) AT,TT,EL1,MM,ZT,ETI,(PT(J),J=1,MM)
      M = min(ns,mm)
      CALL EPTR(EL,EL1,I,*2)
      IF ( I .GT. 0 .AND. IND(I) .EQ. -1) THEN
         ATM(I) = AT
         TRM(I) = TT
         MAX(I) = M
         ZZ(I)  = ZT
         C = 1.d0
         IF ( Z .NE. ZT ) C = Z/ZT
!
!  *****  SCALE RESULTS IF DATA IS FOR AN ATOM WITH A DIFFERENT Z
!
         e(I,I) = C*C*ETI
         DO 11 J = 1,M
            P(J,I) = C*PT(J)
11       CONTINUE
!
!  *****  SET REMAINING VALUES IN THE RANGE = 0.
!
         IF ( M .EQ. ns ) GO TO 12
         M = M +1
         DO 13  J=M,ns
13       P(J,I) = 0.d0
12       IND(I) = -2
      ENDIF
      GO TO 2
!
!  *****  SET PARAMTERS FOR ELECTRONS AND INITIALIZE FUNCTIONS
!
5     continue
      DO 9 I = 1,NWF
      IF (IND(I)) 7,8,9
!
!  ***** WAVE FUNCTIONS NOT FOUND IN THE INPUT DATA, SET IND = 0
!
7     IF ( IND(I) .EQ. -2 ) GO TO 4
      IND(I) = 0
      WRITE(iscw,27) EL(I)
27    FORMAT(8X,'WAVE FUNCTIONS NOT FOUND FOR ',A3)
!
!  *****  DETERMINE ESTIMATES OF THE WAVE FUNCTIONS BY THE SCREENED
!  *****  HYDROGENIC APPROXIMATION
!
8     continue
      CALL BHWF(N(I),L(I),Z-s(i),nt,kx,ks,nint,gr,grw,bsp,bs,w21,w22,
     :          P(1,I))
      M = NS-1
30    IF ( DABS(P(M,I)) .LT. 1.D-15 ) then
        P(M,I) = 0.d0
        M = M-1
        GO TO 30
      END IF
!31    MAX(I) = M+1
 31    MAX(I) = ns
       print *, ' nl expansion', el(i)
       print '(6f12.8)', (p(ii,i),ii=1,ns)
      e(i,i) = ((Z-s(i))/n(i))**2
       print *, 'After orthogonalizationn'
!
!  *****  ORTHOGONALIZE TO INNER FUNCTIONS
!
4     IM = I - 1
      m = max(i)
      DO 6 II =1,IM
	if (e(i,ii) .ne. 0.d0) then
          PN = QUADR(I,II,0)
          print *, 'Overlap between ',i,ii,pn,e(i,ii)
          IF ( DABS(PN) .GT. 1.D-10 ) THEN
            M = MAX0(m,MAX(II))
            DO 25 J = 1,M
 25           P(J,I) =P(J,I) - PN*P(J,II)
          END IF
	end if
6     CONTINUE
      pn = 1.d0/sqrt(quadr(i,i,0))
      if (p(4,i) .lt. 0.d0) pn = -pn
      do 16 j = 1,m
	p(j,i) = pn*p(j,i)
16     continue
      print *, ' expansion for ', n(i),l(i)
      print '(6f12.8)' ,(p(j,i),j=1,ns)
9     CONTINUE
!
!     .. improve estimates obtained from screened hydrogenics
!
      call improve(ind)
!
      WRITE(3,14)
14    FORMAT(/// 8X,18HINITIAL ESTIMATES  //10X,2HNL,
     1   4X,5HSIGMA,6X,5HE(NL),4X,9HFUNCTIONS//)
!
      DO 15 I = 1,NWF
      K = IND(I) + 2
      IF ( IND(I) .EQ. -2 ) THEN
           TITLE = ' SCALED '//ATM(I)//TRM(I)
        ELSE IF (IND(I) .EQ. 0) THEN
           TITLE = ' SCREENED HYDROGENIC'
        ELSE
           TITLE = ' UNCHANGED'
      END IF
17    WRITE(3,19) EL(I),S(I),E(I,I),TITLE
19    FORMAT(9X,A3,F9.2,F11.3,3X,A24)
15    CONTINUE
      ec = 0.d0
      IF (iuf .ne. 0) close(unit=iuf)
!      print *, ' Check orthogonlity'
!      a11 = quadr(1,1,0)
!      a12 = quadr(1,2,0)
!      a22 = quadr(2,2,0)
!      print *, a11, a12, a22
!      print *, '1s', '2s'
!      print '(6f12.8)', (p(ii,1),ii=1,ns)
!      print '(6f12.8)', (p(ii,2),ii=1,ns)
      END
