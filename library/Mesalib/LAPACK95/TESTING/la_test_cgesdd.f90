SUBROUTINE LA_TEST_CGESDD(JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, &
&   RWORK, IWORK, INFO)
!
!  -- LAPACK95 interface driver routine (version 1.1) --
!     UNI-C, Denmark;
!     May 31, 1999
!
!  .. Use Statements ..
      USE LA_PRECISION, ONLY: WP => SP
      USE F95_LAPACK, ONLY: LA_GESDD
!  .. Implicit Statement ..
      IMPLICIT NONE
!  .. Scalar Arguments ..
      INTEGER, INTENT(IN) :: M, N, LDA, LDU, LDVT, LWORK
      INTEGER, INTENT(INOUT) :: INFO
      CHARACTER*1, INTENT(IN) :: JOBZ
!  .. Array Arguments ..
      COMPLEX(WP), INTENT(INOUT) :: A(1:LDA,1:N)
      COMPLEX(WP), INTENT(OUT):: WORK(1:LWORK)
      REAL(WP), INTENT(OUT) :: S(1: MIN(M,N))
      COMPLEX(WP), INTENT(OUT) :: U(1: LDU, 1:N), VT(1: LDVT, 1:N)
      INTEGER :: IWORK(1: 8*MIN(M,N))
      REAL(WP), INTENT (OUT) :: RWORK(1: 5*min(M,N)*min(M,N) + 7*min(M,N))
!  .. Parameters ..
      CHARACTER(LEN=8),  PARAMETER :: SRNAME = 'LA_GESDD'
      CHARACTER(LEN=14), PARAMETER :: SRNAMT = 'LA_TEST_CGESDD'
!  .. Common blocks ..
      INTEGER :: INFOTC
      COMMON /LINFO95/ INFOTC
!  .. Local Scalars ..
      INTEGER :: I, J, IA1, IA2, IS, IU1, IU2, IVT1, IVT2, ISTAT
      CHARACTER*1 IJOBZ, JOB
!  .. Local Arrays ..
      LOGICAL, SAVE :: CTEST = .TRUE., ETEST = .TRUE.
      LOGICAL LSAME
      COMPLEX(WP), POINTER :: W1(:,:), W2(:,:)
      INTEGER MN
!  .. Executable Statements ..
      IA1 = M; IA2 = N; IS = MIN(M,N); IU1 = M; IU2 = M
      IVT1 = N; IVT2 = N; IJOBZ=JOBZ; MN = MIN(M,N)
      I = INFO / 100; J = INFO - I*100
      IF (LSAME(JOBZ, 'S')) THEN
        IU2 = min(M,N)
        IVT2 = min(M,N)
      END IF
      SELECT CASE(I)
      CASE (2)
      IS = IS - 1
      JOB = 'N'
      CALL LA_GESDD( A(1:IA1,1:IA2), S(1:IS), U(1:IU1, 1:IU2), &
     &  VT(1:IVT1, 1:IVT2), JOB=JOB, INFO=INFO)
      CALL LA_AUX_AA01( I, CTEST, ETEST, SRNAMT )
      RETURN
      CASE (3)
      IU1 = IU1 - 1
      JOB = 'N'
      CALL LA_GESDD( A(1:IA1,1:IA2), S(1:IS), U(1:IU1, 1:IU2), &
     &  VT(1:IVT1, 1:IVT2), JOB=JOB, INFO=INFO)
      CALL LA_AUX_AA01( I, CTEST, ETEST, SRNAMT )
      RETURN 
      CASE (4)
      IVT2 = IVT2 - 1
      JOB = 'N'
      CALL LA_GESDD( A(1:IA1,1:IA2), S(1:IS), U(1:IU1, 1:IU2), &
     &  VT(1:IVT1, 1:IVT2), JOB=JOB, INFO=INFO)
      CALL LA_AUX_AA01( I, CTEST, ETEST, SRNAMT )
      RETURN
      CASE (6)
      JOB = 'T'
      CALL LA_GESDD( A(1:IA1,1:IA2), S(1:IS), U(1:IU1, 1:IU2), &
     &  VT(1:IVT1, 1:IVT2), JOB=JOB, INFO=INFO)
      CALL LA_AUX_AA01( I, CTEST, ETEST, SRNAMT )
      RETURN 
      CASE(:-1,1, 5, 7:)
      CALL UESTOP(SRNAMT)
      END SELECT


   SELECT CASE (JOBZ)
     CASE ('A')
       JOB = 'N'
       CALL LA_GESDD( A(1:IA1,1:IA2), S(1:IS), U(1:IU1, 1:IU2), &
     &   VT(1:IVT1, 1:IVT2), JOB=JOB, INFO=INFO)
     CASE ('S')
       JOB = 'N'
       ALLOCATE (W1(M,M), W2(N,N), STAT=ISTAT)
       IF (ISTAT == 0) THEN
         CALL LA_GESDD( A(1:IA1,1:IA2), S(1:IS), W1, &
     &     W2, JOB=JOB, INFO=INFO)
         U(1:M, 1:MN) = W1(1:M, 1:MN)
         VT(1:MN, 1:N) = W2(1:MN, 1:N)
         DEALLOCATE (W2, W1)
       ELSE
         INFO = -111
       ENDIF
       
     CASE ('O')
       IF (M .GE. N) THEN
         JOB = 'U'
         CALL LA_GESDD( A(1:IA1,1:IA2), S(1:IS),  &
     &     VT=VT(1:IVT1, 1:IVT2), JOB=JOB, INFO=INFO)
       ELSE
         JOB = 'V'
         CALL LA_GESDD( A(1:IA1,1:IA2), S(1:IS), U(1:IU1, 1:IU2), &
     &     JOB=JOB, INFO=INFO)
       ENDIF
     CASE ('N')
       JOB = 'N'
       CALL LA_GESDD( A(1:IA1,1:IA2), S(1:IS),  &
&        JOB=JOB, INFO=INFO)
   END SELECT
   
   CALL LA_AUX_AA01( I, CTEST, ETEST, SRNAMT )
 END SUBROUTINE LA_TEST_CGESDD
      
