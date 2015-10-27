      program coulomb  
      implicit none
      integer      MODE,i,KFN,IEXP,IFAIL,LMIN,LSTOP,nmax,k
      double precision RHOMAT,X,XK,RMATCH,ETA,Z,rstep,ra,rmax,energy,fac
      double precision FC(51),GC(51),FCP(51), GCP(51)


!!      MODE = 3
      MODE = 1
      KFN  = 0


      LMIN  = 0
      LSTOP = 0

      rmax = 80.0d0
      ra   = 0.05d0
      rstep = 0.02d0
      nmax  = (rmax-ra)/rstep

      energy  = 1.d0
      XK  =  dsqrt(2.d0*energy)
      Z = 1.d0
      ETA = -Z/XK


         fac = 1.d0
      do i = 0, nmax
         x = ra + dble(i)*rstep

      rhomat = xk*x
      CALL COULFG (RHOMAT,ETA,LMIN,LSTOP,FC,GC,FCP,GCP,MODE,KFN,IEXP,
     +             IFAIL)

      if(IFAIL /= 0) write(*,*) 'Warning message !!!'
      if(i == 0.and.FC(LMIN+1) < 0.d0)  fac = -1.d0 

         FC(LMIN+1)  = fac* FC(LMIN+1)
         FCP(LMIN+1) = fac* FCP(LMIN+1)
         GC(LMIN+1)  = fac* GC(LMIN+1)
         GCP(LMIN+1) = fac* GCP(LMIN+1)

      write(61,400)X,FC(LMIN+1),FCP(LMIN+1),GC(LMIN+1),GCP(LMIN+1)

      end do


100    format(e22.15,3x,e22.15,3x,e22.15)
200    format(e14.8,2x,e12.6,2x,e12.6)
300    format(1x,d22.15,4x,d22.15,4x,d22.15)
400    format(1x,e12.5,4x,e22.15,4x,e22.15, 4x,e22.15,4x,e22.15)
      STOP
      end 

*=======================================================================
*=======================================================================
*                                                                      *
*  C O U L F G  -  P A C K A G E   (FROM  A.R. B A R N E T T)          *
*                                                                      *
*  THIS PACKAGE IS USED TO CALCULATE REGULAR AND IRREGULAR             *
*  COULOMB - AND BESSEL - FUNCTIONS                                    *
*                                                                      *
      SUBROUTINE COULFG(XX,ETA1,LMIN,LMAX,FC,GC,FCP,GCP,
     *                  MODE1,KFN,IEXP,IFAIL)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  REVISED COULOMB WAVEFUNCTION PROGRAM USING STEED'S METHOD           C
C                                                                      C
C  A. R. BARNETT           MANCHESTER  MARCH   1981                    C
C                                                                      C
C  ORIGINAL PROGRAM 'RCWFN'      IN    CPC  8 (1974) 377-395           C
C                 + 'RCWFF'      IN    CPC 11 (1976) 141-142           C
C  FULL DESCRIPTION OF ALGORITHM IN    CPC 21 (1981) 297-314           C
C  THIS VERSION WRITTEN UP       IN    CPC 27 (1982) 147-166           C
C                                                                      C
C                                                                      C
C                                                                      C
C                                                                      C
C  THIS CODE IS A MODIFIED VERSION OF THAT PUBLISHED BY BARNETT IN     C
C  CPC 27 (1982) 147-166. IT HAS BEEN RE-CODED IN FORTRAN 77 AND       C
C  THE ACCURACY IS DETERMINED BY A MACHINE DEPENDENT PARAMETER         C
C  ( SEE BELOW UNDER 'ACCURACY' ).                                     C
C                                                                      C
C                                                                      C
C  COULFG RETURNS F,G,F',G', FOR REAL XX.GT.0,REAL ETA1 (INCLUDING 0), C
C   AND INTEGER LAMBDA.GE.0 FOR A RANGE OF LAMBDA VALUES:              C
C   LMIN TO LMAX.                                                      C
C   STARTING ARRAY ELEMENT IS M1 = LMIN+1                              C
C                                                                      C
C  IF 'MODE' = 1  GET F,G,F',G'                                        C
C            = 2  GET F,G                                              C
C            = 3  GET F                                                C
C  IF 'KFN'  = 0  REAL       COULOMB FUNCTIONS ARE RETURNED            C
C            = 1  SPHERICAL   BESSEL      "      "     "               C
C            = 2  CYLINDRICAL BESSEL      "      "     "               C
C                                                                      C
C  THE USE OF 'MODE' AND 'KFN' IS INDEPENDENT                          C
C                                                                      C
C  ACCURACY                                                            C
C  ========                                                            C
C                                                                      C
C                                                                      C
C  THE ACCURACY IS DETERMINED BY THE VARIABLE: ACCUR                   C
C  ACCUR IS SET TO: MAX(U,1.0D-16), WHERE U IS A MACHINE DEPENDENT     C
C  QUANTITY DETERMINED BY THE SUBROUTINE MACHIN. U IS A MEASURE OF     C
C  THE MACHINE ACCURACY.                                               C
C  THE USER MUST CALL MACHIN BEFORE THE FIRST CALL TO COULFG.          C
C  IF THE USER'S MACHINE ALLOWS MORE PRECISION THAN 1.0D-16 AND IF     C
C  A PRECISION BETTER THAN 1.0D-16 IS REQUIRED, THEN ALTER THE VALUE   C
C  OF 'TM16' IN THE PARAMETER STATEMENT BELOW.                         C
C  IN THE OSCILLATING REGION X.GE.XTURN, WHERE                         C
C  XTURN = ETA1+SQRT( ETA1**2+ LMIN*( LMIN+1) ), SOLUTIONS ARE         C
C  OBTAINED TO AN ACCURACY ACCUR. HOWEVER IF X IS SUFFICIENTLY         C
C  SMALLER THAN XTURN, SO THAT G.GT.1.0D6, THEN SOLUTIONS ARE          C
C  OBTAINED USING A JWKB APPROXIMATION AND THE RESULTS WILL BE MUCH    C
C  LESS ACCURATE ( IN GENERAL THE JWKB APPROXIMATION PROVIDES RESULTS  C
C  TO BETTER THAN 1% ). IF THE JWKB APPROXIMATION IS USED, A WARNING   C
C  MESSAGE IS PRINTED OUT FOR THE USER'S INFORMATION.                  C
C                                                                      C
C                                                                      C
C   OVERFLOW/UNDERFLOW                                                 C
C   ==================                                                 C
C                                                                      C
C   TO AVOID UNDERFLOW/OVERFLOW WHEN THE JWKB APPROXIMATION IS USED    C
C   THE USER MUST SET THE PARAMETER 'IUO' IN SUBROUTINE JWKB.          C
C                                                                      C
C                                                                      C
C   IEXP ON OUTPUT                                                     C
C   ==============                                                     C
C                                                                      C
C   IF IEXP .GT.1 ON OUTPUT, THEN SCALED RESULTS EXIST IN THE ARRAYS   C
C   FC,GC,FCP AND GCP. THE TRUE SOLUTIONS ARE FC*10**(-IEXP);          C
C   FCP*10**(-IEXP); GC*10**(IEXP); AND GCP*10**(IEXP).                C
C                                                                      C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      LOGICAL      ETANE0,XLTURN
      PARAMETER ( IREAD = 5, IWRITE = 6 )
*
*     IREAD AND IWRITE SPECIFY THE INPUT AND OUTPUT UNITS
*
      PARAMETER ( KDIM2 =51 )
*
*     KDIM2 IS THE ARRAY STORAGE DIMENSION FOR THE ARRAYS
*          FC,GC,FCP,GCP
*          IT MUST BE AT LEAST LMAX+1
*
      PARAMETER ( ZERO  = 0.0D0, TM30   = 1.0D-30,  TM16   = 1.0D-16,
     *            HALF  = 0.5D0, ONE    = 1.0D0  ,  TWO    = 2.0D0,
     *            TEN2  = 1.0D2, TEN4   = 1.0D4  ,  TEN6   = 1.0D6,
     *            ABORT = 2.0D+04, ABORT2 = 4.0D4 )
 
      PARAMETER ( RT2DPI =
     *            0.79788456080286535587989211986876373D0)
*
*     RT2DPI = SQRT(TWO/PI)
*
      DIMENSION    FC(KDIM2),GC(KDIM2),FCP(KDIM2),GCP(KDIM2)
      COMMON  / STEED /  PACCQ,NFP,NPQ,M1
*
*     COMMON BLOCK IS FOR INFORMATION ONLY.  NOT REQUIRED IN CODE
*
 
 
      COMMON  / SMALL /  FOURU,TWOU,U
*
*     COULFG HAS CALLS TO: ABS,ANINT,DBLE,MAX,MIN,SIGN,NINT,SQRT
*
*
*     CHECK DIMENSIONS
*
      IF ( KDIM2 .LT. (LMAX+1) ) THEN
        WRITE(IWRITE,1007)
        WRITE(IWRITE,1006) LMAX,KDIM2
        STOP
      ENDIF
      IFAIL = 0
*
*     SET VALUES OF ACCURACY VARIABLES
*
      ACCUR = MAX(U,TM16)
      ACC = ACCUR
      ACC4 = ACC*TEN4
      ACCH = SQRT(ACC)
*
*    TEST RANGE OF XX, EXIT IF.LE.DSQRT(ACCUR) OR IF NEGATIVE
*
      IF ( XX.LE.ACCH ) THEN
        IFAIL = -1
        write(iwrite,*) 'This is IFAIL = -1 :'
        WRITE(IWRITE,1007)
        WRITE(IWRITE,1001) XX,ACCH
        RETURN
      ENDIF
*
*     CHECK THAT LMIN AND LMAX ARE SUCH THAT:
*     XLM.GT.-ONE   AND   LMAX.GE.LMIN
*
      XLM = ANINT(DBLE(LMIN))
      IF ( XLM .LE. -ONE .OR. LMAX .LT. LMIN ) THEN
        IFAIL = -2
        write(iwrite,*) 'This is IFAIL = -2 :'
        WRITE(IWRITE,1007)
        WRITE(IWRITE,1002) LMAX,LMIN,XLM
        RETURN
      ENDIF
      IF ( KFN.EQ.2 ) THEN
        XLM = XLM-HALF
      ENDIF
*
*     DETERMINE LXTRA = THE NUMBER OF ADDITIONAL LAMBDA VALUES
*     TO BE COMPUTED.
*
      LXTRA = LMAX-LMIN
      IF ( MODE1.EQ.1 ) THEN
        MODE = 1
      ELSE
        MODE = MODE1
      ENDIF
      IEXP = 1
      NPQ = 0
      GJWKB = ZERO
      PACCQ = ONE
      IF ( KFN.NE.0 ) THEN
        ETA = ZERO
        ETANE0 = .FALSE.
      ELSE
        ETA = ETA1
        ETANE0 = .TRUE.
      ENDIF
      X = XX
*
*     DETERMINE WHETHER X IS .LT. THE TURNING POINT VALUE
*     IF IT IS: SET XLTURN = .TRUE.
*     IF NOT  : SET XLTURN = .FALSE.
*
      IF ( X*(X - TWO*ETA) .LT. XLM*XLM + XLM ) THEN
        XLTURN = .TRUE.
      ELSE
        XLTURN = .FALSE.
      ENDIF

!!!!      write(*,*) 'XLTURN =', XLTURN


      E2MM1 = ETA*ETA + XLM*XLM + XLM
      XLL = XLM + DBLE(LXTRA)
*
*       XLL  ISMAX LAMBDA VALUE, OR 0.5 SMALLER FOR J,Y BESSELS
*
*         DETERMINE STARTING ARRAY ELEMENT (M1) FROM LMIN
*
      M1 = LMIN+1
      L1 = M1 + LXTRA
*
*    EVALUATE CF1  = F    =  FPRIME(XL,ETA,X)/F(XL,ETA,X)
*
      XI = ONE/X
      FCL = ONE
      PK = XLL + ONE
      PX = PK + ABORT
1     CONTINUE
      EK = ETA / PK
      F = (EK + PK*XI)*FCL + (FCL - ONE)*XI
      PK1 = PK + ONE
*
*   ENSURE THAT B1 .NE. ZERO FOR NEGATIVE ETA: FIXUP IS EXACT.
*
      IF ( ABS(ETA*X + PK*PK1) .GT. ACC ) THEN
        GOTO 2
      ENDIF
      FCL = (ONE + EK*EK)/(ONE + (ETA/PK1)**2)
      PK = TWO + PK
      GOTO 1
2     D = ONE/((PK + PK1)*(XI + EK/PK1))
      DF = -FCL*(ONE + EK*EK)*D
      IF ( FCL .NE. ONE ) THEN
        FCL = -ONE
      ENDIF
      IF ( D .LT. ZERO) THEN
        FCL = -FCL
      ENDIF
      F = F + DF
*
*   BEGIN CF1 LOOP ON PK = K = LAMBDA + 1
*
      P = ONE
3     CONTINUE
      PK = PK1
      PK1 = PK1 + ONE
      EK = ETA / PK
      TK = (PK + PK1)*(XI + EK/PK1)
      D = TK - D*(ONE + EK*EK)

!!!          WRITE(IWRITE,1003) ABORT,F,DF,PK,PX,ACC
!!!          STOP

      IF ( ABS(D) .LE. ACCH ) THEN
        WRITE (IWRITE,1000) D,DF,ACCH,PK,EK,ETA,X
        P = P + ONE
        IF( P .GT. TWO ) THEN
          IFAIL = 1
          write(iwrite,*) 'This is IFAIL = 1 :'
          WRITE(IWRITE,1007)
          WRITE(IWRITE,1003) ABORT,F,DF,PK,PX,ACC
          RETURN
        ENDIF
      ENDIF

      D = ONE/D
      IF ( D .LT. ZERO ) THEN
        FCL = -FCL
      ENDIF
      DF = DF*(D*TK - ONE)
      F = F + DF
      IF ( PK .GT. PX ) THEN
        IFAIL = 1
        write(iwrite,*) 'This is IFAIL = 1 :'
        WRITE(IWRITE,1007)
        WRITE(IWRITE,1003) ABORT,F,DF,PK,PX,ACC
        RETURN
      ENDIF
      IF ( ABS(DF) .LT. ABS(F)*ACC ) THEN
        GOTO 4
      ENDIF
      GOTO 3
4     CONTINUE
      NFP = NINT(PK - XLL - ONE)
      IF ( LXTRA .NE. 0 ) THEN
*
*   DOWNWARD RECURRENCE TO LAMBDA = XLM. ARRAY GC,IF PRESENT,STORES RL
*
        FCL = FCL*TM30
        FPL = FCL*F
        IF ( MODE .EQ. 1 ) THEN
          FCP(L1) = FPL
        ENDIF
        FC (L1) = FCL
        XL = XLL
        RL = ONE
        EL = ZERO
        DO 5, LP = 1,LXTRA
          IF ( ETANE0 ) THEN
            EL = ETA/XL
            RL = SQRT(ONE+EL*EL)
          ENDIF
          SL = EL + XL*XI
          L = L1 - LP
          FCL1 = (FCL *SL + FPL)/RL
          FPL = FCL1*SL - FCL *RL
          FCL = FCL1
          FC(L) = FCL
          IF ( MODE .EQ. 1 ) THEN
            FCP(L) = FPL
          ENDIF
          IF ( MODE .NE. 3 .AND. ETANE0 ) THEN
            GC(L+1) = RL
          ENDIF
          XL = XL - ONE
5       CONTINUE
        IF ( FCL .EQ. ZERO ) THEN
          FCL = ACC
        ENDIF
        F = FPL/FCL
      ENDIF
*
*    NOW WE HAVE REACHED LAMBDA = XLMIN = XLM
*    EVALUATE CF2 = P + I.Q  AGAIN USING STEED'S ALGORITHM
*    SEE TEXT FOR COMPACT COMPLEX CODE FOR SP CDC OR NON-ANSI IBM
*
      IF ( XLTURN ) THEN
        CALL JWKB(X,ETA,MAX(XLM,ZERO),FJWKB,GJWKB,IEXP)
      ENDIF
      IF( IEXP .GT. 1 .OR. GJWKB .GT. TEN6 ) THEN
*
*     AT THIS POINT  G(XLM) .GT. 10**6 OR IEXP .GT. IUO & XLTURN = .TRUE
*
        WRITE(IWRITE,1008)
        W = FJWKB
        GAM = GJWKB*W
        P = F
        Q = ONE
      ELSE
        XLTURN = .FALSE.
        PK = ZERO
        WI = ETA + ETA
        P = ZERO
        Q = ONE - ETA*XI
        AR = -E2MM1
        AI = ETA
        BR = TWO*(X - ETA)
        BI = TWO
        DR = BR/(BR*BR + BI*BI)
        DI = -BI/(BR*BR + BI*BI)
        DP = -XI*(AR*DI + AI*DR)
        DQ = XI*(AR*DR - AI*DI)
6       CONTINUE
        P = P + DP
        Q = Q + DQ
        PK = PK + TWO
        AR = AR + PK
        AI = AI + WI
        BI = BI + TWO
        D = AR*DR - AI*DI + BR
        DI = AI*DR + AR*DI + BI
        C = ONE/(D*D + DI*DI)
        DR = C*D
        DI = -C*DI
        A = BR*DR - BI*DI - ONE
        B = BI*DR + BR*DI
        C = DP*A - DQ*B
        DQ = DP*B + DQ*A
        DP = C
        IF ( PK .GT. ABORT2 ) THEN
          IFAIL = 2
          write(iwrite,*) 'This is IFAIL = 2 :'
          WRITE(IWRITE,1007)
          WRITE(IWRITE,1004) ABORT,P,Q,DP,DQ,ACC
          RETURN
        ENDIF
        IF ( ABS(DP)+ABS(DQ).LT.(ABS(P)+ABS(Q))*ACC ) THEN
          GOTO 7
        ENDIF
        GOTO 6
7       CONTINUE
        NPQ = NINT(PK/TWO)
        PACCQ = HALF*ACC/MIN(ABS(Q),ONE)
        IF ( ABS(P).GT.ABS(Q) ) THEN
          PACCQ = PACCQ*ABS(P)
        ENDIF
        GAM = (F - P)/Q
        IF ( Q .LE. ACC4*ABS(P) ) THEN
          IFAIL = 3
          write(iwrite,*) 'This is IFAIL = 3 :'
          WRITE(IWRITE,1007)
          WRITE(IWRITE,1005) P,Q,ACC,LXTRA,M1
          RETURN
        ENDIF
        W = ONE/SQRT((F - P)*GAM + Q)
      ENDIF
*
*    NORMALISE FOR SPHERICAL OR CYLINDRICAL BESSEL FUNCTIONS
*
      IF ( KFN.EQ.0 ) THEN
        ALPHA = ZERO
        BETA = ONE
      ELSE
     +  IF ( KFN.EQ.1 ) THEN
          ALPHA = XI
          BETA = XI
        ELSE
          ALPHA = XI*HALF
          BETA = SQRT(XI)*RT2DPI
      ENDIF
      FCM = SIGN(W,FCL)*BETA
      FC(M1) = FCM
      IF ( MODE.LT.3 ) THEN
        IF ( XLTURN ) THEN
          GCL = GJWKB*BETA
        ELSE
          GCL = FCM*GAM
        ENDIF
        IF ( KFN.NE.0 ) THEN
          GCL = -GCL
        ENDIF
        GC(M1) = GCL
        GPL = GCL*(P - Q/GAM) - ALPHA*GCL
        IF ( MODE.EQ.1 ) THEN
          GCP(M1) = GPL
          FCP(M1) = FCM*(F - ALPHA)
        ENDIF
      ENDIF
      IF ( LXTRA .NE. 0 ) THEN
*
*     UPWARD RECURRENCE FROM GC(M1),GCP(M1)  STORED VALUE IS RL
*     RENORMALISE FC,FCP AT EACH LAMBDA AND CORRECT REGULAR DERIVATIVE
*        XL   = XLM HERE  AND RL = ONE , EL = ZERO FOR BESSELS
*
        W = BETA*W/ABS(FCL)
        MAXL = L1 - 1
        DO 8, L = M1,MAXL
          IF ( MODE .LT. 3 ) THEN
            XL = XL + ONE
            IF ( ETANE0 ) THEN
              EL = ETA/XL
              RL = GC(L+1)
            ENDIF
            SL = EL + XL*XI
            GCL1 = ((SL - ALPHA)*GCL - GPL)/RL
            GPL = RL*GCL - (SL + ALPHA)*GCL1
            GCL = GCL1
            GC(L+1) = GCL1
            IF ( MODE .EQ. 1 ) THEN
              GCP(L+1) = GPL
              FCP(L+1) = W*(FCP(L+1) - ALPHA*FC(L+1))
            ENDIF
          ENDIF
          FC(L+1) = W* FC(L+1)
8       CONTINUE
      ENDIF
      RETURN
*
*     FORMAT STATEMENTS:
*
1000  FORMAT (/' CF1 ACCURACY LOSS: D,DF,ACCH,K,ETA/K,ETA,X = ',7D
     *         9.2/)
1001  FORMAT(' FOR XX = ',D12.3,' TRY SMALL-X  SOLUTIONS',
     *' OR X NEGATIVE'/ ,' SQUARE ROOT ACCURACY PARAMETER = ',D12.3/)
1002  FORMAT(/' PROBLEM WITH INPUT ORDER VALUES:LMAX,LMIN,XLM = ',
     *2I6,D15.6/)
1003  FORMAT(' CF1 HAS FAILED TO CONVERGE AFTER ',F10.0,' ITERATIONS',/
     *' F,DF,PK,PX,ACCUR =  ',5D12.3//)
1004  FORMAT(' CF2 HAS FAILED TO CONVERGE AFTER ',F7.0,' ITERATIONS',/
     *' P,Q,DP,DQ,ACCUR =  ',4D17.7,D12.3//)
1005  FORMAT(' FINAL Q.LE.DABS(P)*ACC*10**4 , P,Q,ACC = ',3D12.3,4X,
     *' LXTRA,M1 = ',2I5 /)
1006  FORMAT(' ARRAY DIMENSIONS SHOULD BE INCREASED TO INT(LMAX+1).'/
     +  ' COULFG WAS CALLED WITH LMAX = ',I6/
     +  ' PRESENT DIMENSIONS ARE KDIM2 = ',I6)
1007  FORMAT(//' COULFG FAILURE. '//)
1008  FORMAT(//' INFORMATION MESSAGE.'/
     *' SOLUTIONS WERE OBTAINED USING THE JWKB APPROXIMATION '/
     *' THEY MAY BE LESS ACCURATE THAN DESIRED'//)
      END
* ==============================================================
      SUBROUTINE JWKB(XX,ETA1,XL,FJWKB,GJWKB,IEXP)
      DOUBLE PRECISION          XX,ETA1,XL,FJWKB,GJWKB,DZERO
*
*     COMPUTES JWKB APPROXIMATIONS TO COULOMB FUNCTIONS   FOR XL.GE. 0
*     AS MODIFIED BY BIEDENHARN ET AL. PHYS REV 97 (1955) 542-554
*
*     CALLS ATAN2,EXP,INT,LOG,MAX,REAL,SQRT
*
 
 
 
 
      PARAMETER ( IUO = 70 )
*
*    IUO MUST BE SET BY THE USER. ITS PURPOSE IS TO PREVENT
*    UNDERFLOW/OVERFLOW AND SO IT IS MACHINE DEPENDENT.
*    IUO SHOULD HAVE AS VALUE THE EXPONENT OF THE LARGEST
*    REAL NUMBER WHICH THE MACHINE IS CAPABLE OF HOLDING
*    WITHOUT OVERFLOW MINUS 5.
*
*
      PARAMETER ( ZERO  = 0.0E0, HALF  = 0.5E0, ONE  = 1.0E0,
     *            SIX   = 6.0E0, TEN   = 1.0E1, RL35  = 3.5E1,
     *            ALOGE = 0.4342945E0 )
 
      PARAMETER ( DZERO = 0.0D0 )
 
 
      X     = XX
      ETA   = ETA1
      GH2   = X*(ETA + ETA - X)
      XLL1  = MAX(XL*XL + XL,DZERO)
      IF ( GH2 + XLL1 .GT. ZERO ) THEN
        HLL = XLL1 + SIX/RL35
        HL = SQRT(HLL)
        SL = ETA/HL + HL/X
        RL2 = ONE + ETA*ETA/HLL
        GH = SQRT(GH2 + HLL)/X
        PHI = X*GH - HALF*( HL*LOG((GH + SL)**2/RL2) - LOG(GH) )
        IF ( ETA.NE.ZERO ) THEN
          PHI = PHI - ETA*ATAN2(X*GH,X - ETA)
        ENDIF
        PHI10 = -PHI*ALOGE
        IEXP = INT(PHI10)
        IF ( IEXP.GT.IUO ) THEN
          GJWKB = TEN**(PHI10 - REAL(IEXP))
        ELSE
          GJWKB = EXP(-PHI)
          IEXP = 0
        ENDIF
        FJWKB = HALF/(GH*GJWKB)
      ENDIF
      RETURN
      END
*====================================================================
      SUBROUTINE MACHIN
*
*
*     U IS THE SMALLEST POSITIVE NUMBER SUCH THAT
*     (1.+U) .GT. 1.
*     U IS COMPUTED APPROXIMATELY AS A POWER OF 1./2.
*
*     THIS ROUTINE IS COMPLETELY EXPLAINED AND DOCUMENTED
*     IN THE TEXT:
*     " COMPUTER SOLUTION OF ORDINARY DIFFERENTIAL EQUATIONS:
*       THE INITIAL VALUE PROBLEM " BY L.F.SHAMPINE AND M.K.GORDON
*
*
*
*     THE ROUTINE HAS BEEN RE-CODED IN FORTRAN 77
*
*
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER ( HALF=0.5D0 , ONE=1.0D0 , TWO=2.0D0 )
      COMMON / SMALL / FOURU,TWOU,U
      HALFU = HALF
1     T = ONE+HALFU
      IF ( T.LE.ONE ) THEN
        U = TWO*HALFU
        TWOU = TWO*U
        FOURU = TWO*TWOU
      ELSE
        HALFU = HALF*HALFU
        GOTO 1
      ENDIF
      RETURN
      END


