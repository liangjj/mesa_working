      SUBROUTINE ZDRVPP( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX,
     $                   A, AFAC, ASAV, B, BSAV, X, XACT, S, WORK,
     $                   RWORK, NOUT )
*
*  -- LAPACK test routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            TSTERR
      INTEGER            NMAX, NN, NOUT, NRHS
      DOUBLE PRECISION   THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      INTEGER            NVAL( * )
      DOUBLE PRECISION   RWORK( * ), S( * )
      COMPLEX*16         A( * ), AFAC( * ), ASAV( * ), B( * ),
     $                   BSAV( * ), WORK( * ), X( * ), XACT( * )
*     ..
*
*  Purpose
*  =======
*
*  ZDRVPP tests the driver routines ZPPSV and -SVX.
*
*  Arguments
*  =========
*
*  DOTYPE  (input) LOGICAL array, dimension (NTYPES)
*          The matrix types to be used for testing.  Matrices of type j
*          (for 1 <= j <= NTYPES) are used for testing if DOTYPE(j) =
*          .TRUE.; if DOTYPE(j) = .FALSE., then type j is not used.
*
*  NN      (input) INTEGER
*          The number of values of N contained in the vector NVAL.
*
*  NVAL    (input) INTEGER array, dimension (NN)
*          The values of the matrix dimension N.
*
*  NRHS    (input) INTEGER
*          The number of right hand side vectors to be generated for
*          each linear system.
*
*  THRESH  (input) DOUBLE PRECISION
*          The threshold value for the test ratios.  A result is
*          included in the output file if RESULT >= THRESH.  To have
*          every test ratio printed, use THRESH = 0.
*
*  TSTERR  (input) LOGICAL
*          Flag that indicates whether error exits are to be tested.
*
*  NMAX    (input) INTEGER
*          The maximum value permitted for N, used in dimensioning the
*          work arrays.
*
*  A       (workspace) COMPLEX*16 array, dimension (NMAX*(NMAX+1)/2)
*
*  AFAC    (workspace) COMPLEX*16 array, dimension (NMAX*(NMAX+1)/2)
*
*  ASAV    (workspace) COMPLEX*16 array, dimension (NMAX*(NMAX+1)/2)
*
*  B       (workspace) COMPLEX*16 array, dimension (NMAX*NRHS)
*
*  BSAV    (workspace) COMPLEX*16 array, dimension (NMAX*NRHS)
*
*  X       (workspace) COMPLEX*16 array, dimension (NMAX*NRHS)
*
*  XACT    (workspace) COMPLEX*16 array, dimension (NMAX*NRHS)
*
*  S       (workspace) DOUBLE PRECISION array, dimension (NMAX)
*
*  WORK    (workspace) COMPLEX*16 array, dimension
*                      (NMAX*max(3,NRHS))
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (NMAX+2*NRHS)
*
*  NOUT    (input) INTEGER
*          The unit number for output.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      INTEGER            NTYPES
      PARAMETER          ( NTYPES = 9 )
      INTEGER            NTESTS
      PARAMETER          ( NTESTS = 6 )
*     ..
*     .. Local Scalars ..
      LOGICAL            EQUIL, NOFACT, PREFAC, ZEROT
      CHARACTER          DIST, EQUED, FACT, PACKIT, TYPE, UPLO, XTYPE
      CHARACTER*3        PATH
      INTEGER            I, IEQUED, IFACT, IMAT, IN, INFO, IOFF, IUPLO,
     $                   IZERO, K, K1, KL, KU, LDA, MODE, N, NERRS,
     $                   NFACT, NFAIL, NIMAT, NPP, NRUN, NT
      DOUBLE PRECISION   AINVNM, AMAX, ANORM, CNDNUM, RCOND, RCONDC,
     $                   ROLDC, SCOND
*     ..
*     .. Local Arrays ..
      CHARACTER          EQUEDS( 2 ), FACTS( 3 ), PACKS( 2 ), UPLOS( 2 )
      INTEGER            ISEED( 4 ), ISEEDY( 4 )
      DOUBLE PRECISION   RESULT( NTESTS )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DGET06, ZLANSP
      EXTERNAL           LSAME, DGET06, ZLANSP
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALADHD, ALAERH, ALASVM, ZCOPY, ZERRVX, ZGET04,
     $                   ZLACPY, ZLAQSP, ZLARHS, ZLASET, ZLATB4, ZLATMS,
     $                   ZPPEQU, ZPPSV, ZPPSVX, ZPPT01, ZPPT02, ZPPT05,
     $                   ZPPTRF, ZPPTRI
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      CHARACTER*6        SRNAMT
      INTEGER            INFOT, NUNIT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCMPLX, MAX
*     ..
*     .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS / 'U', 'L' / , FACTS / 'F', 'N', 'E' / ,
     $                   PACKS / 'C', 'R' / , EQUEDS / 'N', 'Y' /
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      PATH( 1: 1 ) = 'Zomplex precision'
      PATH( 2: 3 ) = 'PP'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
*
*     Test the error exits
*
      IF( TSTERR )
     $   CALL ZERRVX( PATH, NOUT )
      INFOT = 0
*
*     Do for each value of N in NVAL
*
      DO 150 IN = 1, NN
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         NPP = N*( N+1 ) / 2
         XTYPE = 'N'
         NIMAT = NTYPES
         IF( N.LE.0 )
     $      NIMAT = 1
*
         DO 140 IMAT = 1, NIMAT
*
*           Do the tests only if DOTYPE( IMAT ) is true.
*
            IF( .NOT.DOTYPE( IMAT ) )
     $         GO TO 140
*
*           Skip types 3, 4, or 5 if the matrix size is too small.
*
            ZEROT = IMAT.GE.3 .AND. IMAT.LE.5
            IF( ZEROT .AND. N.LT.IMAT-2 )
     $         GO TO 140
*
*           Do first for UPLO = 'U', then for UPLO = 'L'
*
            DO 130 IUPLO = 1, 2
               UPLO = UPLOS( IUPLO )
               PACKIT = PACKS( IUPLO )
*
*              Set up parameters with ZLATB4 and generate a test matrix
*              with ZLATMS.
*
               CALL ZLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE,
     $                      CNDNUM, DIST )
               RCONDC = ONE / CNDNUM
*
               SRNAMT = 'ZLATMS'
               CALL ZLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE,
     $                      CNDNUM, ANORM, KL, KU, PACKIT, A, LDA, WORK,
     $                      INFO )
*
*              Check error code from ZLATMS.
*
               IF( INFO.NE.0 ) THEN
                  CALL ALAERH( PATH, 'ZLATMS', INFO, 0, UPLO, N, N, -1,
     $                         -1, -1, IMAT, NFAIL, NERRS, NOUT )
                  GO TO 130
               END IF
*
*              For types 3-5, zero one row and column of the matrix to
*              test that INFO is returned correctly.
*
               IF( ZEROT ) THEN
                  IF( IMAT.EQ.3 ) THEN
                     IZERO = 1
                  ELSE IF( IMAT.EQ.4 ) THEN
                     IZERO = N
                  ELSE
                     IZERO = N / 2 + 1
                  END IF
*
*                 Set row and column IZERO of A to 0.
*
                  IF( IUPLO.EQ.1 ) THEN
                     IOFF = ( IZERO-1 )*IZERO / 2
                     DO 20 I = 1, IZERO - 1
                        A( IOFF+I ) = ZERO
   20                CONTINUE
                     IOFF = IOFF + IZERO
                     DO 30 I = IZERO, N
                        A( IOFF ) = ZERO
                        IOFF = IOFF + I
   30                CONTINUE
                  ELSE
                     IOFF = IZERO
                     DO 40 I = 1, IZERO - 1
                        A( IOFF ) = ZERO
                        IOFF = IOFF + N - I
   40                CONTINUE
                     IOFF = IOFF - IZERO
                     DO 50 I = IZERO, N
                        A( IOFF+I ) = ZERO
   50                CONTINUE
                  END IF
               ELSE
                  IZERO = 0
               END IF
*
*              Save a copy of the matrix A in ASAV.
*
               CALL ZCOPY( NPP, A, 1, ASAV, 1 )
*
               DO 120 IEQUED = 1, 2
                  EQUED = EQUEDS( IEQUED )
                  IF( IEQUED.EQ.1 ) THEN
                     NFACT = 3
                  ELSE
                     NFACT = 1
                  END IF
*
                  DO 110 IFACT = 1, NFACT
                     FACT = FACTS( IFACT )
                     PREFAC = LSAME( FACT, 'F' )
                     NOFACT = LSAME( FACT, 'N' )
                     EQUIL = LSAME( FACT, 'E' )
*
                     IF( ZEROT ) THEN
                        IF( PREFAC )
     $                     GO TO 110
                        RCONDC = ZERO
*
                     ELSE IF( .NOT.LSAME( FACT, 'N' ) ) THEN
*
*                       Compute the condition number for comparison with
*                       the value returned by ZPPSVX (FACT = 'N' reuses
*                       the condition number from the previous iteration
*                          with FACT = 'F').
*
                        CALL ZCOPY( NPP, ASAV, 1, AFAC, 1 )
                        IF( EQUIL .OR. IEQUED.GT.1 ) THEN
*
*                          Compute row and column scale factors to
*                          equilibrate the matrix A.
*
                           CALL ZPPEQU( UPLO, N, AFAC, S, SCOND, AMAX,
     $                                  INFO )
                           IF( INFO.EQ.0 .AND. N.GT.0 ) THEN
                              IF( IEQUED.GT.1 )
     $                           SCOND = ZERO
*
*                             Equilibrate the matrix.
*
                              CALL ZLAQSP( UPLO, N, AFAC, S, SCOND,
     $                                     AMAX, EQUED )
                           END IF
                        END IF
*
*                       Save the condition number of the
*                       non-equilibrated system for use in ZGET04.
*
                        IF( EQUIL )
     $                     ROLDC = RCONDC
*
*                       Compute the 1-norm of A.
*
                        ANORM = ZLANSP( '1', UPLO, N, AFAC, RWORK )
*
*                       Factor the matrix A.
*
                        CALL ZPPTRF( UPLO, N, AFAC, INFO )
*
*                       Form the inverse of A.
*
                        CALL ZCOPY( NPP, AFAC, 1, A, 1 )
                        CALL ZPPTRI( UPLO, N, A, INFO )
*
*                       Compute the 1-norm condition number of A.
*
                        AINVNM = ZLANSP( '1', UPLO, N, A, RWORK )
                        IF( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
                           RCONDC = ONE
                        ELSE
                           RCONDC = ( ONE / ANORM ) / AINVNM
              Ba-$!j46%c0M%`-$P&$58a-cBa-$!j46%
c0M%`-$P&-63f06!`18%a0$Be-$!j36%d0M8`-$P"-63f06!`18%a0$Be-$!j36%
e0MN0*6!`16Ba06Bj-$!j0M%e0MN`-$Nf-68f16!`16Ba0$Be-$!j36%d0M8`-$P
"-63f06!`18%a0$Be-$!j33dP-63f06!`18%a0$Be-$!j36%d0M8`-$P"-6-f-6!
`188a-cBa-$!j46%c083`-%%b-6-e4$!`36)a-c9%$58`-%%b-6)e16!`36Ba-M8
j-$""0M%c083`-%%b-6-f-6!`188a-cBa-$!j46%d0M8`-$P"-63f06!`18%0*6%
e0MN`-$Nf-6Bf4$!`16)a0MC%-$!j-M%f0c%`-$K&-6Fh06!`1%%a1$Fj-$!i0M%
i0cN`-$Jf-6Nh4!dP-$!i-M%j0d3`-$Jb-6Ni-6!`0d8a16Ja-$!h46%j0d3`-$J
b-6Nh4$!`1$)a1$Fj-$!i0M%i0cN`-$Jf$58a0cFe-$!i36%f0c%`-$K&-6Bf4$!
`16)a0MC%-$!j-M%e0MN`-$Nf-68f16!`16Ba06Bj-$!j0M%d0M80*6!`18%a0$B
e-$!j36%d0M8`-$P"-63f06!`18%a0$Be-$!j36%f0N3`-$Nb-6Fh06!`1%%a16G
%-$!i-JdP-6Ni-6!`0d8a3MJj-$!h0M&$16%`-$C&-88j16!`0MBa4MP%-$!f-M)
`368`-$9"-M*"36!`068b-d&&$58`-$8a-M4#0M!`0$Nb0N*&-$!d-6)f3c)`-$0
%-MG$0M!`-cNb0N-b-$!c4$)e3N%`-$3e-M0"46!`06%0*6)`368`-$9"-88j16!
`0MBa3cNa-$!f46&"1$8`-$G"-6Nh4$!`1$)a0cFe-$!i36%f0N3`-$Nb-63f03d
P-$!j36%c0M%`-$P&-6)e16!`36B`4M8`-$""4M"'0%-`-%)c-%3d0$!`3N)`4$3
d-$"#3M"&0$J`-%)h$58`4M4$-$"#-c"'0%-`-%)c-%Be-$!`38Ba-68e-$""36%
b06N`-%%f-6-e4$!`36)a-cBa-$!j46%d0M80*6!`18%a0$Be-$!j36%e0MN`-$N
f-6Bf4$!`16)a0MFa-$!i46%h0c8`-$K"-6Jh16!`1$Ba1$Fj-$!i0JdP-6Fh06!
`1%%a0MFa-$!i46%f0c%`-$K&-6Bf4$!`16)a0MC%-$!j-M%e0MN`-$Nf-63f06!
`18%a0$Be$58`-$P"-63f06!`18%a0$Be-$!j36%d0M8`-$P"-6-f-6!`188a-cB
a-$!j46%c0M%`-$P&-6-e4$!`36)0*6%c083`-%%b-6-f-6!`188a-cBa-$!j46%
d0M8`-$P"-63f06!`18%a0$Be-$!j36%e0MN`-$Nf-68f13dP-$!j0M%e0MN`-$N
f-68f16!`16Ba06Bj-$!j0M%e0MN`-$Nf-6Bf4$!`16)a0MC%-$!j-M%f0N3`-$N
b$58a0MC%-$!j-M%e0MN`-$Nf-68f16!`16Ba06Bj-$!j0M%d0M8`-$P"-63f06!
`18%a0$Be-$!j36%d0M80*6!`18%a-cBa-$!j46%c0M%`-$P&-6-f-6!`188a-cB
a-$!j46%d0M8`-$P"-63f06!`18%a0$Be-$!j33dP-63f06!`18%a0$Be-$!j36%
d0M8`-$P"-63f06!`18%a0$Be-$!j36%d0M8`-$P"-63f06!`18%a0$Be$58`-$P
"-68f16!`16Ba06Bj-$!j0M%e0MN`-$Nf-63f06!`18%a0$Be-$!j36%d0M8`-$P
"-63f06!`18%0*6%c0M%`-$P&-6-f-6!`188a-cBa-$!j46%c0M%`-$P&-6-f-6!
`188a-cBa-$!j46%c0M%`-$P&-6-f-3dP-$!j46%c0M%`-$P&-6-f-6!`188a-cB
a-$!j46%d0M8`-$P"-63f06!`18%a0$Be-$!j36%d0M8`-$P"$58a0$Be-$!j36%
d0M8`-$P"-63f06!`18%a06Bj-$!j0M%e0MN`-$Nf-63f06!`18%a0$Be-$!j36%
d0M80*6!`18%a0$Be-$!j36%d0M8`-$P"-63f06!`18%a0$Be-$!j36%d0M8`-$P
"-63f06!`18%a0$Be-$!j33dP-63f06!`18%a0$Be-$!j36%d0M8`-$P"-63f06!
`18%a0$Be-$!j36%d0M8`-$P"-63f06!`18%a0$Be$58`-$P"-63f06!`18%a0$B
e-$!j36%d0M8`-$P"-63f06!`18%a0$Be-$!j36%d0M8`-$P"-63f06!`18%0*6%
d0M8`-$P"-63f06!`18%a0$Be-$!j36%d0M8`-$P"-63f06!`18%a0$Be-$!j36%
d0M8`-$P"-63f03dP-$!j36%d0M8`-$P"-63f06!`18%a0$Be-$!j36%c0M%`-$P
&-6-f-6!`188a-cBa-$!j46%c0M%`-$P&$58a-cBa-$!j46%c0M%`-$P&-6-f-6!
`188a-cBa-$!j46%c0M%`-$P&-6-f-6!`188a-cBa-$!j46%c0M%0*6!`188a0$B
e-$!j36%d0M8`-$P"-63f06!`18%a0$Be-$!j36%d0M8`-$P"-68f16!`16Ba06B
j-$!j0JdP-63f06!`18%a0$Be-$!j36%d0M8`-$P"-63f06!`18%a0$Be-$!j36%
c0M%`-$P&-6-f-6!`188a-cBa$58`-$P&-6-f-6!`188a-cBa-$!j46%c0M%`-$P
&-6-f-6!`188a-cBa-$!j46%c0M%`-$P&-6-f-6!`1880*6%c0M%`-$P&-6-f-6!
`188a0$Be-$!j36%d0M8`-$P"-63f06!`18%a0$Be-$!j36%e0MN`-$Nf-68f13d
P-$!j0M%e0MN`-$Nf-68f16!`16Ba06Bj-$!j0M%d0M8`-$P"-63f06!`18%a0$B
e-$!j36%d0M8`-$P"$58a0$Be-$!j36%d0M8`-$P"-6-f-6!`188a-cBa-$!j46%
c0M%`-$P&-6-e4$!`36)a-c9%-$""-M%b06N0*6!`36Ba-M8j-$""0M%b06N`-%%
f-6-e4$!`36)a-cBa-$!j46%d0M8`-$P"-63f06!`18%a0$Be-$!j33dP-68f16!
`16Ba0MC%-$!j-M%f0c%`-$K&-6Bh-6!`1%8a0cFe-$!i36%i0cN`-$Jf-6Nh4$!
`1$)a16G%$58`-$Jb-6Nh4$!`1$)a16Ja-$!h46%j1$%`-$G&-6Nh4$!`1$)a16G
%-$!i-M%i0cN`-$Jf-6Fh06!`1%%0*6%f0c%`-$K&-6Bf4$!`16)a06Bj-$!j0M%
e0MN`-$Nf-63f06!`18%a0$Be-$!j36%d0M8`-$P"-63f03dP-$!j36%d0M8`-$P
"-6-f-6!`188a-cBa-$!j46%d0M8`-$P"-68f16!`16Ba0MFa-$!i46%i0cN`-$J
f$58a16Ja-$!h46&#1$N`-$Ff-8-j-6!`0N8a4$Ne-$!f36&'183`-$Bb-M""06!
`08%b-N&"-$!e06)c3880*6!`06%b-d)b-$!d4$)e3N%`-$3e-MC#46!`0$%b0N-
b-$!c4$)f3c)`-$0%-M4#0M!`0$Nb-N&"-$!e03dP-8C"-6!`088a4$Ne-$!f36&
$1%3`-$Fb-8%i06!`0d%a16G%-$!i-M%f0c%`-$K&-68f16!`16Ba0$Be$58`-$P
"-6-e4$!`36)a-68e-$""36"'0%-`-%)c-%8d1$!`3MF`3c3`-$"#4M"$-d-`-%-
c-%-d-$!`3NB0*6"&0$J`-%)h-%Bd3c!`3M-`4M4$-$"#-c"'06!`-%&'-6%e06!
`38%a-M8j-$""0M%c083`-%%b-6-f-3dP-$!j46%d0M8`-$P"-63f06!`18%a06B
j-$!j0M%f0c%`-$K&-6Fh06!`1%%a1$Fj-$!i0M%i0cN`-$Jf$58a0cFe-$!i36%
h0c8`-$K"-6Bh-6!`1%8a0MC%-$!j-M%f0N3`-$Nb-68f16!`16Ba06Bj-$!j0M%
d0M80*6!`18%a0$Be-$!j36%d0M8`-$P"-63f06!`18%a-cBa-$!j46%c0M%`-$P
&-6-f-6!`188a-c9%-$""-JdP-6-e4$!`36)a-cBa-$!j46%c0M%`-$P&-63f06!
`18%a0$Be-$!j36%d0M8`-$P"-68f16!`16Ba06Bj$58`-$Nf-68f16!`16Ba06B
j-$!j0M%e0MN`-$Nf-6Bf4$!`16)a0MC%-$!j-M%f0N3`-$Nb-6Bf4$!`16)0*6%
f0N3`-$Nb-6Bf4$!`16)a06Bj-$!j0M%e0MN`-$Nf-63f06!`18%a0$Be-$!j36%
d0M8`-$P"-63f03dP-$!j36%c0M%`-$P&-6-f-6!`188a-cBa-$!j46%c0M%`-$P
&-6-f-6!`188a0$Be-$!j36%d0M8`-$P"$58a0$Be-$!j36%d0M8`-$P"-63f06!
`18%a0$Be-$!j36%d0M8`-$P"-63f06!`18%a0$Be-$!j36%d0M80*6!`18%a06B
j-$!j0M%e0MN`-$Nf-68f16!`16Ba0$Be-$!j36%d0M8`-$P"-63f06!`18%a0$B
e-$!j33dP-6-f-6!`188a-cBa-$!j46%c0M%`-$P&-6-f-6!`188a-cBa-$!j46%
c0M%`-$P&-6-f-6!`188a-cBa$58`-$P&-6-f-6!`188a-cBa-$!j46%c0M%`-$P
&-63f06!`18%a0$Be-$!j36%d0M8`-$P"-63f06!`18%0*6%d0M8`-$P"-63f06!
`18%a0$Be-$!j36%e0MN`-$Nf-63f06!`18%a0$Be-$!j36%d0M8`-$P"-63f03d
P-$!j36%d0M8`-$P"-63f06!`18%a0$Be-$!j36%d0M8`-$P"-63f06!`18%a0$B
e-$!j36%d0M8`-$P"$58a0$Be-$!j36%d0M8`-$P"-63f06!`18%a0$Be-$!j36%
d0M8`-$P"-63f06!`18%a0$Be-$!j36%d0M80*6!`18%a0$Be-$!j36%d0M8`-$P
"-63f06!`18%a0$Be-$!j36%d0M8`-$P"-63f06!`18%a0$Be-$!j33dP-63f06!
`18%a0$Be-$!j36%d0M8`-$P"-63f06!`18%a0$Be-$!j36%d0M8`-$P"-63f06!
`18%a0$Be$58`-$P"-63f06!`18%a0$Be-$!j36%d0M8`-$P"-6-f-6!`188a-cB
a-$!j46%c0M%`-$P&-6-f-6!`1880*6%c0M%`-$P&-6-f-6!`188a-cBa-$!j46%
c0M%`-$P&-6-f-6!`188a-cBa-$!j46%c0M%`-$P&-6-f-3dP-$!j46%d0M8`-$P
"-63f06!`18%a0$Be-$!j36%d0M8`-$P"-63f06!`18%a06Bj-$!j0M%e0MN`-$N
f$58a0$Be-$!j36%d0M8`-$P"-63f06!`18%a0$Be-$!j36%d0M8`-$P"-6-f-6!
`188a-cB            ROLDC, RESULT( 3 ) )
                        END IF
*
*                       Check the error bounds from iterative
*                       refinement.
*
                        CALL ZPPT05( UPLO, N, NRHS, ASAV, B, LDA, X,
     $                               LDA, XACT, LDA, RWORK,
     $                               RWORK( NRHS+1 ), RESULT( 4 ) )
                     ELSE
                        K1 = 6
                     END IF
*
*                    Compare RCOND from ZPPSVX with the computed value
*                    in RCONDC.
*
                     RESULT( 6 ) = DGET06( RCOND, RCONDC )
*
*                    Print information about the tests that did not pass
*                    the threshold.
*
                     DO 90 K = K1, 6
                        IF( RESULT( K ).GE.THRESH ) THEN
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                        CALL ALADHD( NOUT, PATH )
                           IF( PREFAC ) THEN
                              WRITE( NOUT, FMT = 9997 )'ZPPSVX', FACT,
     $                           UPLO, N, EQUED, IMAT, K, RESULT( K )
                           ELSE
                              WRITE( NOUT, FMT = 9998 )'ZPPSVX', FACT,
     $                           UPLO, N, IMAT, K, RESULT( K )
                           END IF
                           NFAIL = NFAIL + 1
                        END IF
   90                CONTINUE
                     NRUN = NRUN + 7 - K1
  100                CONTINUE
  110             CONTINUE
  120          CONTINUE
  130       CONTINUE
  140    CONTINUE
  150 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASVM( PATH, NOUT, NFAIL, NRUN, NERRS )
*
 9999 FORMAT( 1X, A6, ', UPLO=''', A1, ''', N =', I5, ', type ', I1,
     $      ', test(', I1, ')=', G12.5 )
 9998 FORMAT( 1X, A6, ', FACT=''', A1, ''', UPLO=''', A1, ''', N=', I5,
     $      ', type ', I1, ', test(', I1, ')=', G12.5 )
 9997 FORMAT( 1X, A6, ', FACT=''', A1, ''', UPLO=''', A1, ''', N=', I5,
     $      ', EQUED=''', A1, ''', type ', I1, ', test(', I1, ')=',
     $      G12.5 )
      RETURN
*
*     End of ZDRVPP
*
      END
