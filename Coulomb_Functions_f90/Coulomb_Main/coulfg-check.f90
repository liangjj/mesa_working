PROGRAM coulomb
 
! Code converted using TO_F90 by Alan Miller
! Date: 2008-01-11  Time: 14:40:29
 
IMPLICIT NONE
INTEGER :: mode,i,kfn,iexp,ifail,lmin,lstop,nmax,k
DOUBLE PRECISION :: rhomat,x,xk,rmatch,eta,z,rstep,ra,rmax,energy,fac
DOUBLE PRECISION :: fc(51),gc(51),fcp(51), gcp(51)


!!      MODE = 3
mode = 1
kfn  = 0


lmin  = 0
lstop = 0

rmax = 80.0D0
ra   = 0.05D0
rstep = 0.02D0
nmax  = (rmax-ra)/rstep

energy  = 1.d0
xk  =  DSQRT(2.d0*energy)
z = 1.d0
eta = -z/xk


fac = 1.d0
DO i = 0, nmax
  x = ra + DBLE(i)*rstep
  
  rhomat = xk*x
  CALL coulfg (rhomat,eta,lmin,lstop,fc,gc,fcp,gcp,mode,kfn,iexp, ifail)
  
  IF(ifail /= 0) WRITE(*,*) 'Warning message !!!'
  IF(i == 0.AND.fc(lmin+1) < 0.d0)  fac = -1.d0
  
  fc(lmin+1)  = fac* fc(lmin+1)
  fcp(lmin+1) = fac* fcp(lmin+1)
  gc(lmin+1)  = fac* gc(lmin+1)
  gcp(lmin+1) = fac* gcp(lmin+1)
  
  WRITE(61,400)x,fc(lmin+1),fcp(lmin+1),gc(lmin+1),gcp(lmin+1)
  
END DO


100    FORMAT(e22.15,3X,e22.15,3X,e22.15)
200    FORMAT(e14.8,2X,e12.6,2X,e12.6)
300    FORMAT(1X,d22.15,4X,d22.15,4X,d22.15)
400    FORMAT(1X,e12.5,4X,e22.15,4X,e22.15, 4X,e22.15,4X,e22.15)
STOP
END PROGRAM coulomb

!=======================================================================
!=======================================================================
!                                                                      *
!  C O U L F G  -  P A C K A G E   (FROM  A.R. B A R N E T T)          *
!                                                                      *
!  THIS PACKAGE IS USED TO CALCULATE REGULAR AND IRREGULAR             *
!  COULOMB - AND BESSEL - FUNCTIONS                                    *
!                                                                      *

SUBROUTINE coulfg(xx,eta1,lmin,lmax,fc,gc,fcp,gcp, mode1,kfn,iexp,ifail)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                      C
!  REVISED COULOMB WAVEFUNCTION PROGRAM USING STEED'S METHOD           C
!                                                                      C
!  A. R. BARNETT           MANCHESTER  MARCH   1981                    C
!                                                                      C
!  ORIGINAL PROGRAM 'RCWFN'      IN    CPC  8 (1974) 377-395           C
!                 + 'RCWFF'      IN    CPC 11 (1976) 141-142           C
!  FULL DESCRIPTION OF ALGORITHM IN    CPC 21 (1981) 297-314           C
!  THIS VERSION WRITTEN UP       IN    CPC 27 (1982) 147-166           C
!                                                                      C
!                                                                      C
!                                                                      C
!                                                                      C
!  THIS CODE IS A MODIFIED VERSION OF THAT PUBLISHED BY BARNETT IN     C
!  CPC 27 (1982) 147-166. IT HAS BEEN RE-CODED IN FORTRAN 77 AND       C
!  THE ACCURACY IS DETERMINED BY A MACHINE DEPENDENT PARAMETER         C
!  ( SEE BELOW UNDER 'ACCURACY' ).                                     C
!                                                                      C
!                                                                      C
!  COULFG RETURNS F,G,F',G', FOR REAL XX.GT.0,REAL ETA1 (INCLUDING 0), C
!   AND INTEGER LAMBDA.GE.0 FOR A RANGE OF LAMBDA VALUES:              C
!   LMIN TO LMAX.                                                      C
!   STARTING ARRAY ELEMENT IS M1 = LMIN+1                              C
!                                                                      C
!  IF 'MODE' = 1  GET F,G,F',G'                                        C
!            = 2  GET F,G                                              C
!            = 3  GET F                                                C
!  IF 'KFN'  = 0  REAL       COULOMB FUNCTIONS ARE RETURNED            C
!            = 1  SPHERICAL   BESSEL      "      "     "               C
!            = 2  CYLINDRICAL BESSEL      "      "     "               C
!                                                                      C
!  THE USE OF 'MODE' AND 'KFN' IS INDEPENDENT                          C
!                                                                      C
!  ACCURACY                                                            C
!  ========                                                            C
!                                                                      C
!                                                                      C
!  THE ACCURACY IS DETERMINED BY THE VARIABLE: ACCUR                   C
!  ACCUR IS SET TO: MAX(U,1.0D-16), WHERE U IS A MACHINE DEPENDENT     C
!  QUANTITY DETERMINED BY THE SUBROUTINE MACHIN. U IS A MEASURE OF     C
!  THE MACHINE ACCURACY.                                               C
!  THE USER MUST CALL MACHIN BEFORE THE FIRST CALL TO COULFG.          C
!  IF THE USER'S MACHINE ALLOWS MORE PRECISION THAN 1.0D-16 AND IF     C
!  A PRECISION BETTER THAN 1.0D-16 IS REQUIRED, THEN ALTER THE VALUE   C
!  OF 'TM16' IN THE PARAMETER STATEMENT BELOW.                         C
!  IN THE OSCILLATING REGION X.GE.XTURN, WHERE                         C
!  XTURN = ETA1+SQRT( ETA1**2+ LMIN*( LMIN+1) ), SOLUTIONS ARE         C
!  OBTAINED TO AN ACCURACY ACCUR. HOWEVER IF X IS SUFFICIENTLY         C
!  SMALLER THAN XTURN, SO THAT G.GT.1.0D6, THEN SOLUTIONS ARE          C
!  OBTAINED USING A JWKB APPROXIMATION AND THE RESULTS WILL BE MUCH    C
!  LESS ACCURATE ( IN GENERAL THE JWKB APPROXIMATION PROVIDES RESULTS  C
!  TO BETTER THAN 1% ). IF THE JWKB APPROXIMATION IS USED, A WARNING   C
!  MESSAGE IS PRINTED OUT FOR THE USER'S INFORMATION.                  C
!                                                                      C
!                                                                      C
!   OVERFLOW/UNDERFLOW                                                 C
!   ==================                                                 C
!                                                                      C
!   TO AVOID UNDERFLOW/OVERFLOW WHEN THE JWKB APPROXIMATION IS USED    C
!   THE USER MUST SET THE PARAMETER 'IUO' IN SUBROUTINE JWKB.          C
!                                                                      C
!                                                                      C
!   IEXP ON OUTPUT                                                     C
!   ==============                                                     C
!                                                                      C
!   IF IEXP .GT.1 ON OUTPUT, THEN SCALED RESULTS EXIST IN THE ARRAYS   C
!   FC,GC,FCP AND GCP. THE TRUE SOLUTIONS ARE FC*10**(-IEXP);          C
!   FCP*10**(-IEXP); GC*10**(IEXP); AND GCP*10**(IEXP).                C
!                                                                      C
!                                                                      C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

IMPLICIT DOUBLE PRECISION (a-h,o-z)
INTEGER, PARAMETER :: kdim2 =51
DOUBLE PRECISION, INTENT(IN)             :: xx
DOUBLE PRECISION, INTENT(IN)             :: eta1
INTEGER, INTENT(IN)                      :: lmin
INTEGER, INTENT(IN)                      :: lmax
DOUBLE PRECISION, INTENT(OUT)            :: fc(kdim2)
DOUBLE PRECISION, INTENT(OUT)            :: gc(kdim2)
DOUBLE PRECISION, INTENT(OUT)            :: fcp(kdim2)
DOUBLE PRECISION, INTENT(OUT)            :: gcp(kdim2)
INTEGER, INTENT(IN)                      :: mode1
INTEGER, INTENT(IN)                      :: kfn
INTEGER, INTENT(OUT)                     :: iexp
INTEGER, INTENT(OUT)                     :: ifail
SAVE
LOGICAL :: etane0,xlturn
INTEGER, PARAMETER :: iread = 5
INTEGER, PARAMETER :: iwrite = 6

!     IREAD AND IWRITE SPECIFY THE INPUT AND OUTPUT UNITS



!     KDIM2 IS THE ARRAY STORAGE DIMENSION FOR THE ARRAYS
!          FC,GC,FCP,GCP
!          IT MUST BE AT LEAST LMAX+1

DOUBLE PRECISION, PARAMETER :: zero  = 0.0D0
DOUBLE PRECISION, PARAMETER :: tm30   = 1.0D-30
DOUBLE PRECISION, PARAMETER :: tm16   = 1.0D-16
DOUBLE PRECISION, PARAMETER :: half  = 0.5D0
DOUBLE PRECISION, PARAMETER :: one    = 1.0D0
DOUBLE PRECISION, PARAMETER :: two    = 2.0D0
DOUBLE PRECISION, PARAMETER :: ten2  = 1.0D2
DOUBLE PRECISION, PARAMETER :: ten4   = 1.0D4
DOUBLE PRECISION, PARAMETER :: ten6   = 1.0D6
DOUBLE PRECISION, PARAMETER :: abort = 2.0D+04
DOUBLE PRECISION, PARAMETER :: abort2 = 4.0D4

DOUBLE PRECISION, PARAMETER :: rt2dpi =0.79788456080286535587989211986876373d0

!     RT2DPI = SQRT(TWO/PI)


COMMON  / steed /  paccq,nfp,npq,m1

!     COMMON BLOCK IS FOR INFORMATION ONLY.  NOT REQUIRED IN CODE



COMMON  / small /  fouru,twou,u

!     COULFG HAS CALLS TO: ABS,ANINT,DBLE,MAX,MIN,SIGN,NINT,SQRT


!     CHECK DIMENSIONS

IF ( kdim2 < (lmax+1) ) THEN
  WRITE(iwrite,1007)
  WRITE(iwrite,1006) lmax,kdim2
  STOP
END IF
ifail = 0

!     SET VALUES OF ACCURACY VARIABLES

accur = MAX(u,tm16)
acc = accur
acc4 = acc*ten4
acch = SQRT(acc)

!    TEST RANGE OF XX, EXIT IF.LE.DSQRT(ACCUR) OR IF NEGATIVE

IF ( xx <= acch ) THEN
  ifail = -1
  WRITE(iwrite,*) 'This is IFAIL = -1 :'
  WRITE(iwrite,1007)
  WRITE(iwrite,1001) xx,acch
  RETURN
END IF

!     CHECK THAT LMIN AND LMAX ARE SUCH THAT:
!     XLM.GT.-ONE   AND   LMAX.GE.LMIN

xlm = ANINT(DBLE(lmin))
IF ( xlm <= -one .OR. lmax < lmin ) THEN
  ifail = -2
  WRITE(iwrite,*) 'This is IFAIL = -2 :'
  WRITE(iwrite,1007)
  WRITE(iwrite,1002) lmax,lmin,xlm
  RETURN
END IF
IF ( kfn == 2 ) THEN
  xlm = xlm-half
END IF

!     DETERMINE LXTRA = THE NUMBER OF ADDITIONAL LAMBDA VALUES
!     TO BE COMPUTED.

lxtra = lmax-lmin
IF ( mode1 == 1 ) THEN
  mode = 1
ELSE
  mode = mode1
END IF
iexp = 1
npq = 0
gjwkb = zero
paccq = one
IF ( kfn /= 0 ) THEN
  eta = zero
  etane0 = .false.
ELSE
  eta = eta1
  etane0 = .true.
END IF
x = xx

!     DETERMINE WHETHER X IS .LT. THE TURNING POINT VALUE
!     IF IT IS: SET XLTURN = .TRUE.
!     IF NOT  : SET XLTURN = .FALSE.

IF ( x*(x - two*eta) < xlm*xlm + xlm ) THEN
  xlturn = .true.
ELSE
  xlturn = .false.
END IF

!!!!      write(*,*) 'XLTURN =', XLTURN


e2mm1 = eta*eta + xlm*xlm + xlm
xll = xlm + DBLE(lxtra)

!       XLL  ISMAX LAMBDA VALUE, OR 0.5 SMALLER FOR J,Y BESSELS

!         DETERMINE STARTING ARRAY ELEMENT (M1) FROM LMIN

m1 = lmin+1
l1 = m1 + lxtra

!    EVALUATE CF1  = F    =  FPRIME(XL,ETA,X)/F(XL,ETA,X)

xi = one/x
fcl = one
pk = xll + one
px = pk + abort
1     CONTINUE
ek = eta / pk
f = (ek + pk*xi)*fcl + (fcl - one)*xi
pk1 = pk + one

!   ENSURE THAT B1 .NE. ZERO FOR NEGATIVE ETA: FIXUP IS EXACT.

IF ( ABS(eta*x + pk*pk1) > acc ) THEN
  GO TO 2
END IF
fcl = (one + ek*ek)/(one + (eta/pk1)**2)
pk = two + pk
GO TO 1
2     d = one/((pk + pk1)*(xi + ek/pk1))
df = -fcl*(one + ek*ek)*d
IF ( fcl /= one ) THEN
  fcl = -one
END IF
IF ( d < zero) THEN
  fcl = -fcl
END IF
f = f + df

!   BEGIN CF1 LOOP ON PK = K = LAMBDA + 1

p = one
3     CONTINUE
pk = pk1
pk1 = pk1 + one
ek = eta / pk
tk = (pk + pk1)*(xi + ek/pk1)
d = tk - d*(one + ek*ek)

!!!          WRITE(IWRITE,1003) ABORT,F,DF,PK,PX,ACC
!!!          STOP

IF ( ABS(d) <= acch ) THEN
  WRITE (iwrite,1000) d,df,acch,pk,ek,eta,x
  p = p + one
  IF( p > two ) THEN
    ifail = 1
    WRITE(iwrite,*) 'This is IFAIL = 1 :'
    WRITE(iwrite,1007)
    WRITE(iwrite,1003) abort,f,df,pk,px,acc
    RETURN
  END IF
END IF

d = one/d
IF ( d < zero ) THEN
  fcl = -fcl
END IF
df = df*(d*tk - one)
f = f + df
IF ( pk > px ) THEN
  ifail = 1
  WRITE(iwrite,*) 'This is IFAIL = 1 :'
  WRITE(iwrite,1007)
  WRITE(iwrite,1003) abort,f,df,pk,px,acc
  RETURN
END IF
IF ( ABS(df) < ABS(f)*acc ) THEN
  GO TO 4
END IF
GO TO 3
4     CONTINUE
nfp = nint(pk - xll - one)
IF ( lxtra /= 0 ) THEN
  
!   DOWNWARD RECURRENCE TO LAMBDA = XLM. ARRAY GC,IF PRESENT,STORES RL
  
  fcl = fcl*tm30
  fpl = fcl*f
  IF ( mode == 1 ) THEN
    fcp(l1) = fpl
  END IF
  fc (l1) = fcl
  xl = xll
  rl = one
  el = zero
  DO   lp = 1,lxtra
    IF ( etane0 ) THEN
      el = eta/xl
      rl = SQRT(one+el*el)
    END IF
    sl = el + xl*xi
    l = l1 - lp
    fcl1 = (fcl *sl + fpl)/rl
    fpl = fcl1*sl - fcl *rl
    fcl = fcl1
    fc(l) = fcl
    IF ( mode == 1 ) THEN
      fcp(l) = fpl
    END IF
    IF ( mode /= 3 .AND. etane0 ) THEN
      gc(l+1) = rl
    END IF
    xl = xl - one
  END DO
  IF ( fcl == zero ) THEN
    fcl = acc
  END IF
  f = fpl/fcl
END IF

!    NOW WE HAVE REACHED LAMBDA = XLMIN = XLM
!    EVALUATE CF2 = P + I.Q  AGAIN USING STEED'S ALGORITHM
!    SEE TEXT FOR COMPACT COMPLEX CODE FOR SP CDC OR NON-ANSI IBM

IF ( xlturn ) THEN
  CALL jwkb(x,eta,MAX(xlm,zero),fjwkb,gjwkb,iexp)
END IF
IF( iexp > 1 .OR. gjwkb > ten6 ) THEN
  
!     AT THIS POINT  G(XLM) .GT. 10**6 OR IEXP .GT. IUO & XLTURN = .TRUE
  
  WRITE(iwrite,1008)
  w = fjwkb
  gam = gjwkb*w
  p = f
  q = one
ELSE
  xlturn = .false.
  pk = zero
  wi = eta + eta
  p = zero
  q = one - eta*xi
  ar = -e2mm1
  ai = eta
  br = two*(x - eta)
  bi = two
  dr = br/(br*br + bi*bi)
  di = -bi/(br*br + bi*bi)
  dp = -xi*(ar*di + ai*dr)
  dq = xi*(ar*dr - ai*di)
  6       CONTINUE
  p = p + dp
  q = q + dq
  pk = pk + two
  ar = ar + pk
  ai = ai + wi
  bi = bi + two
  d = ar*dr - ai*di + br
  di = ai*dr + ar*di + bi
  c = one/(d*d + di*di)
  dr = c*d
  di = -c*di
  a = br*dr - bi*di - one
  b = bi*dr + br*di
  c = dp*a - dq*b
  dq = dp*b + dq*a
  dp = c
  IF ( pk > abort2 ) THEN
    ifail = 2
    WRITE(iwrite,*) 'This is IFAIL = 2 :'
    WRITE(iwrite,1007)
    WRITE(iwrite,1004) abort,p,q,dp,dq,acc
    RETURN
  END IF
  IF ( ABS(dp)+ABS(dq) < (ABS(p)+ABS(q))*acc ) THEN
    GO TO 7
  END IF
  GO TO 6
  7       CONTINUE
  npq = nint(pk/two)
  paccq = half*acc/MIN(ABS(q),one)
  IF ( ABS(p) > ABS(q) ) THEN
    paccq = paccq*ABS(p)
  END IF
  gam = (f - p)/q
  IF ( q <= acc4*ABS(p) ) THEN
    ifail = 3
    WRITE(iwrite,*) 'This is IFAIL = 3 :'
    WRITE(iwrite,1007)
    WRITE(iwrite,1005) p,q,acc,lxtra,m1
    RETURN
  END IF
  w = one/SQRT((f - p)*gam + q)
END IF

!    NORMALISE FOR SPHERICAL OR CYLINDRICAL BESSEL FUNCTIONS

IF ( kfn == 0 ) THEN
  alpha = zero
  beta = one
ELSE  &
  IF ( kfn == 1 ) THEN
    alpha = xi
    beta = xi
  ELSE
    alpha = xi*half
    beta = SQRT(xi)*rt2dpi
  END IF
  fcm = SIGN(w,fcl)*beta
  fc(m1) = fcm
  IF ( mode < 3 ) THEN
    IF ( xlturn ) THEN
      gcl = gjwkb*beta
    ELSE
      gcl = fcm*gam
    END IF
    IF ( kfn /= 0 ) THEN
      gcl = -gcl
    END IF
    gc(m1) = gcl
    gpl = gcl*(p - q/gam) - alpha*gcl
    IF ( mode == 1 ) THEN
      gcp(m1) = gpl
      fcp(m1) = fcm*(f - alpha)
    END IF
  END IF
  IF ( lxtra /= 0 ) THEN
    
!     UPWARD RECURRENCE FROM GC(M1),GCP(M1)  STORED VALUE IS RL
!     RENORMALISE FC,FCP AT EACH LAMBDA AND CORRECT REGULAR DERIVATIVE
!        XL   = XLM HERE  AND RL = ONE , EL = ZERO FOR BESSELS
    
    w = beta*w/ABS(fcl)
    maxl = l1 - 1
    DO   l = m1,maxl
      IF ( mode < 3 ) THEN
        xl = xl + one
        IF ( etane0 ) THEN
          el = eta/xl
          rl = gc(l+1)
        END IF
        sl = el + xl*xi
        gcl1 = ((sl - alpha)*gcl - gpl)/rl
        gpl = rl*gcl - (sl + alpha)*gcl1
        gcl = gcl1
        gc(l+1) = gcl1
        IF ( mode == 1 ) THEN
          gcp(l+1) = gpl
          fcp(l+1) = w*(fcp(l+1) - alpha*fc(l+1))
        END IF
      END IF
      fc(l+1) = w* fc(l+1)
    END DO
  END IF
  RETURN
  
!     FORMAT STATEMENTS:
  
  1000  FORMAT (/' CF1 ACCURACY LOSS: D,DF,ACCH,K,ETA/K,ETA,X = ',7D 9.2/)
  1001  FORMAT(' FOR XX = ',d12.3,' TRY SMALL-X  SOLUTIONS',  &
      ' OR X NEGATIVE'/ ,' SQUARE ROOT ACCURACY PARAMETER = ',d12.3/)
  1002  FORMAT(/' PROBLEM WITH INPUT ORDER VALUES:LMAX,LMIN,XLM = ',  &
      2I6,d15.6/)
  1003  FORMAT(' CF1 HAS FAILED TO CONVERGE AFTER ',f10.0,' ITERATIONS',/  &
      ' F,DF,PK,PX,ACCUR =  ',5D12.3//)
  1004  FORMAT(' CF2 HAS FAILED TO CONVERGE AFTER ',f7.0,' ITERATIONS',/  &
      ' P,Q,DP,DQ,ACCUR =  ',4D17.7,d12.3//)
  1005  FORMAT(' FINAL Q.LE.DABS(P)*ACC*10**4 , P,Q,ACC = ',3D12.3,4X,  &
      ' LXTRA,M1 = ',2I5 /)
  1006  FORMAT(' ARRAY DIMENSIONS SHOULD BE INCREASED TO INT(LMAX+1).'/  &
      ' COULFG WAS CALLED WITH LMAX = ',i6/ ' PRESENT DIMENSIONS ARE KDIM2 = ',i6)
  1007  FORMAT(//' COULFG FAILURE. '//)
  1008  FORMAT(//' INFORMATION MESSAGE.'/  &
      ' SOLUTIONS WERE OBTAINED USING THE JWKB APPROXIMATION '/  &
      ' THEY MAY BE LESS ACCURATE THAN DESIRED'//)
END SUBROUTINE coulfg
! ==============================================================

SUBROUTINE jwkb(xx,eta1,xl,fjwkb,gjwkb,iexp)

DOUBLE PRECISION, INTENT(IN)             :: xx
DOUBLE PRECISION, INTENT(IN)             :: eta1
DOUBLE PRECISION, INTENT(IN)             :: xl
DOUBLE PRECISION, INTENT(OUT)            :: fjwkb
DOUBLE PRECISION, INTENT(OUT)            :: gjwkb
INTEGER, INTENT(OUT)                     :: iexp


!     COMPUTES JWKB APPROXIMATIONS TO COULOMB FUNCTIONS   FOR XL.GE. 0
!     AS MODIFIED BY BIEDENHARN ET AL. PHYS REV 97 (1955) 542-554

!     CALLS ATAN2,EXP,INT,LOG,MAX,REAL,SQRT





INTEGER, PARAMETER :: iuo = 70

!    IUO MUST BE SET BY THE USER. ITS PURPOSE IS TO PREVENT
!    UNDERFLOW/OVERFLOW AND SO IT IS MACHINE DEPENDENT.
!    IUO SHOULD HAVE AS VALUE THE EXPONENT OF THE LARGEST
!    REAL NUMBER WHICH THE MACHINE IS CAPABLE OF HOLDING
!    WITHOUT OVERFLOW MINUS 5.


REAL, PARAMETER :: zero  = 0.0E0
REAL, PARAMETER :: half  = 0.5E0
REAL, PARAMETER :: one  = 1.0E0
REAL, PARAMETER :: six   = 6.0E0
REAL, PARAMETER :: ten   = 1.0E1
REAL, PARAMETER :: rl35  = 3.5E1
REAL, PARAMETER :: aloge = 0.4342945E0

DOUBLE PRECISION, PARAMETER :: dzero = 0.0D0


x     = xx
eta   = eta1
gh2   = x*(eta + eta - x)
xll1  = MAX(xl*xl + xl,dzero)
IF ( gh2 + xll1 > zero ) THEN
  hll = xll1 + six/rl35
  hl = SQRT(hll)
  sl = eta/hl + hl/x
  rl2 = one + eta*eta/hll
  gh = SQRT(gh2 + hll)/x
  phi = x*gh - half*( hl*LOG((gh + sl)**2/rl2) - LOG(gh) )
  IF ( eta /= zero ) THEN
    phi = phi - eta*ATAN2(x*gh,x - eta)
  END IF
  phi10 = -phi*aloge
  iexp = INT(phi10)
  IF ( iexp > iuo ) THEN
    gjwkb = ten**(phi10 - REAL(iexp))
  ELSE
    gjwkb = EXP(-phi)
    iexp = 0
  END IF
  fjwkb = half/(gh*gjwkb)
END IF
RETURN
END SUBROUTINE jwkb
!====================================================================

SUBROUTINE machin


!     U IS THE SMALLEST POSITIVE NUMBER SUCH THAT
!     (1.+U) .GT. 1.
!     U IS COMPUTED APPROXIMATELY AS A POWER OF 1./2.

!     THIS ROUTINE IS COMPLETELY EXPLAINED AND DOCUMENTED
!     IN THE TEXT:
!     " COMPUTER SOLUTION OF ORDINARY DIFFERENTIAL EQUATIONS:
!       THE INITIAL VALUE PROBLEM " BY L.F.SHAMPINE AND M.K.GORDON



!     THE ROUTINE HAS BEEN RE-CODED IN FORTRAN 77



IMPLICIT DOUBLE PRECISION (a-h,o-z)
SAVE
PARAMETER ( half=0.5D0 , one=1.0D0 , two=2.0D0 )
COMMON / small / fouru,twou,u
halfu = half
1     t = one+halfu
IF ( t <= one ) THEN
  u = two*halfu
  twou = two*u
  fouru = two*twou
ELSE
  halfu = half*halfu
  GO TO 1
END IF
RETURN
END SUBROUTINE machin


