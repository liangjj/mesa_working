*     ==================================================================
*    		W A V E F N
*     ==================================================================
*
*       This routine initializes radial functions by the procedure
*   indicated by IND(I).
*
*         Value of IND(I)     Method
*         ---------------     ------
*             -1           Functions read from unit IU2
*              0           Screened hydrogenic functions with ZZ=Z-S(I)
*              1           Functions in memory left unchanged
*                                                  0
*   The set of functions are then orthogonalized, Y (i, i;r) and the
*   diagonal energy parameters computed, when necessary.
*
*
      SUBROUTINE WAVEFN(ind)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER (nt=100,kx=8,nwd=10)
*
      INTEGER ouf
      COMMON /INOUT/ iuf,ouf,iscw
*
      CHARACTER config*50,el*3,atom*6,term*6
      COMMON /LABEL/config,el(nwd),atom,term
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
*
      COMMON /SPLINE/gr(nt,kx),grm(nt,kx),grw(nt,kx),
     :        bsp(nt,kx,kx),bspd(nt,kx,kx,2),db1(nt,kx),db2(nt,kx),
     :        sb(nt,kx),r1(nt,kx),rm1(nt,kx),rm2(nt,kx),
     :        h0(nt,kx), rkm(nt,nt,kx,kx), lh, krk
*
      LOGICAL fail,omit,rel,all,trace,print,strong
      COMMON /TEST/fail,omit,rel,all,trace
*
      DIMENSION ZZ(NWD),PT(NT),ind(nwd),w3(nt),bs(kx,nt),
     :          w21(nt,kx),w22(nt,kx)
      CHARACTER*3 EL1
      CHARACTER*6 AT,TT,ATM(NWD),TRM(NWD)
      CHARACTER*24 TITLE
      LOGICAL LD
*
      nint=ns-ks+1
      call facsb(nt,kx,ks,ns,sb,bs)
*
*  ***** READ THE WAVEFUNCTIONS
*
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
*
*  *****  SCALE RESULTS IF DATA IS FOR AN ATOM WITH A DIFFERENT Z
*
         e(I,I) = C*C*ETI
         DO 11 J = 1,M
            P(J,I) = C*PT(J)
11       CONTINUE
*
*  *****  SET REMAINING VALUES IN THE RANGE = 0.
*
         IF ( M .EQ. ns ) GO TO 12
         M = M +1
         DO 13  J=M,ns
13       P(J,I) = 0.d0
12       IND(I) = -2
      ENDIF
      GO TO 2
*
*  *****  SET PARAMTERS FOR ELECTRONS AND INITIALIZE FUNCTIONS
*
5     continue
      DO 9 I = 1,NWF
      IF (IND(I)) 7,8,9
*
*  ***** WAVE FUNCTIONS NOT FOUND IN THE INPUT DATA, SET IND = 0
*
7     IF ( IND(I) .EQ. -2 ) GO TO 4
      IND(I) = 0
      WRITE(iscw,27) EL(I)
27    FORMAT(8X,'WAVE FUNCTIONS NOT FOUND FOR ',A3)
*
*  *****  DETERMINE ESTIMATES OF THE WAVE FUNCTIONS BY THE SCREENED
*  *****  HYDROGENIC APPROXIMATION
*
8     continue
      CALL BHWF(N(I),L(I),Z-s(i),nt,kx,ks,nint,gr,grw,bsp,bs,w21,w22,
     :          P(1,I))
      M = NS-1
30    IF ( DABS(P(M,I)) .LT. 1.D-15 ) then
        P(M,I) = 0.d0
        M = M-1
        GO TO 30
      END IF
*31    MAX(I) = M+1
 31    MAX(I) = ns
       print *, ' nl expansion', el(i)
       print '(6f12.8)', (p(ii,i),ii=1,ns)
      e(i,i) = ((Z-s(i))/n(i))**2
       print *, 'After orthogonalizationn'
*
*  *****  ORTHOGONALIZE TO INNER FUNCTIONS
*
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
      WRITE(3,14)
14    FORMAT(/// 8X,18HINITIAL ESTIMATES  //10X,2HNL,
     1   4X,5HSIGMA,6X,5HE(NL),4X,9HFUNCTIONS//)
*
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
      END
