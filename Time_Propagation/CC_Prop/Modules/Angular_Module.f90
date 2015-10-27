!***********************************************************************
                           MODULE Angular_Module
                           USE dvr_shared
                           USE dvr_global
                           USE dvrprop_global
                           IMPLICIT NONE
!***********************************************************************
!***********************************************************************
                           Contains
!***********************************************************************
!***********************************************************************
!deck Three_J_Coef
!***begin prologue     Three_J_Coef
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           spherical harmonic
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            Compute 3-J coefficient
!***                   
!***description        See below
!***references
!***routines called    
!***end prologue       Three_J_Coef
!***********************************************************************
  FUNCTION Three_J_Coef(L1,M1,L2,M2,L,M,clebsch)
  IMPLICIT NONE
  INTEGER                                 :: L1, L2, M1, M2, L ,M
  INTEGER, DIMENSION(6)                   :: idum
  REAL*8                                  :: Three_J_Fac, Three_J_Coef
  REAL*8, External                        :: Z_3j
  LOGICAL                                 :: clebsch
!
!        enter as 2*value + 1 to routine z_3j
!
  idum(1) = L1 + L1 + 1
  idum(2) = M1 + M1 + 1
  idum(3) = L2 + L2 + 1
  idum(4) = M2 + M2 + 1
  idum(5) = L + L + 1
  idum(6) = M + M + 1
  IF ( clebsch ) THEN
       Three_J_Fac=(-1)**( ( idum(1) - idum(3) + idum(6) - 1 ) / 2 ) * sqrt( DBLE(idum(5)) )
       Three_J_Coef = Three_J_Fac * Z_3j( idum(1),idum(2),idum(3),idum(4),idum(5), - idum(6) + 2 )
  ELSE
       Three_J_Coef = Z_3j( idum(1),idum(2),idum(3),idum(4),idum(5), idum(6) )
  END IF
!
!***********************************************************************
  END FUNCTION Three_J_Coef
!***********************************************************************
!deck Six_J_Coef
!***begin prologue     Six_J_Coef
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           spherical harmonic
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            Compute 6-J coefficient
!***                   
!***description        See below
!***references
!***routines called    
!***end prologue       Six_J_Coef
!***********************************************************************
  FUNCTION Six_J_Coef(L1,L2,L3,L4,L5,L6)
  IMPLICIT NONE
  INTEGER                                 :: L1, L2, L3, L4, L5 ,L6
  INTEGER, DIMENSION(6)                   :: idum
  REAL*8                                  :: Six_J_Coef
  REAL*8, External                        :: Z_6j
! 
!    enter as 2*value + 1 to routine z_6j
!
  idum(1) = L1 + L1 + 1
  idum(2) = L2 + L2 + 1
  idum(3) = L3 + L3 + 1
  idum(4) = L4 + L4 + 1
  idum(5) = L5 + L5 + 1
  idum(6) = L6 + L6 + 1
  Six_J_Coef = Z_6j( idum(1),idum(2),idum(3),idum(4),idum(5),idum(6))
!
!***********************************************************************
  END FUNCTION Six_J_Coef
!***********************************************************************
!***********************************************************************
  FUNCTION Z_3j (j1,m1,j2,m2,j3,m3) 
  Implicit NONE 
  Real*8                              :: Z_3j 
  Integer, intent(in)                 :: j1,m1,j2,m2,j3,m3 
  Integer                             :: i,i_max,k,kk,m,iz,iz_min,iz_max 
  Real*8                              :: x,y,z 
  Integer, DIMENSION(3)               :: a,b 
  Integer, DIMENSION(16)              :: J 
!-------------------------------------------------------------------- 
! 
!     determines the value of the 3j-symbols without direct using of 
!     factorials. The following expression for the 3j-symbols is used: 
!         (A.P.JUCYS, A.A.BANDZAITIS, 1977) 
! 
!     3j{j1,m1,j2,m2,j3,m3} = delta(m1+m2,m3) * (2j3+1)^1/2 * {j1,j2,j3} * 
!       sqrt[ (j1+m1)!*(j1-m1)!*(j2+m2)!*(j2-m2)!*(j3+m3)!*(j3-m3)! ] 
!                         SUM(z) {   (-1)^z  / 
!          [ z! *  (j1+j2-j3-z)! * (j1-m1-z)! * (j2-m2-z)! * 
!                  (j3-j2+m1+z)! * (j3-j1-m2+z)! ] } 
! 
!     where {a,b,c}=sqrt[ (a+b-c)! * (a-b+c)! * (b-a+c)! / (a+b+c+1)! ] 
! 
!     If we introduce the auxiliary values a(i) and b(i) 
!     (see below the text of program) then 
! 
!     3j =         (-1) ^ Sum[a(i)] 
!          sqrt{ Pr[ (b(j)-a(i))! ] / [ Sum (b(j)-a(i))+1 ] } 
!                  Sum(z) { (-1)^z  / 
!          [  Pr[ (z-a(i))! ]  * Pr[ (b(j)-z)! ]   ] } 
! 
!     (below the moments are used in (2J+1)-representation) 
! 
!-------------------------------------------------------------------- 
  Z_3j=0.d0 
  IF(M1+M2+M3-3 /= 0) RETURN    ! check of conservation rules 
  J(1)= J1+J2-J3-1 
  J(2)= J1-J2+J3-1 
  J(3)= J2-J1+J3-1 
  J(4)= J1+M1-2 
  J(5)= J1-M1 
  J(6)= J2-M2 
  J(7)= J2+M2-2 
  J(8)= J3+M3-2 
  J(9)= J3-M3 
  Do I=1,9 
     IF(J(i) < 0.or.mod(J(i),2) == 1) THEN
        RETURN 
     END IF 
  End do 
  a(1) = 0                         ! auxiliary values 
  a(2) = (j2-j3-m1+1)/2 
  a(3) = (j1-j3+m2-1)/2 
  b(1) = (j1+j2-j3-1)/2 
  b(2) = (j1-m1)/2 
  b(3) = (j2+m2-2)/2 
  IZ_min=MAX0(a(1),a(2),a(3))      ! limits of the sum 
  IZ_max=MIN0(b(1),b(2),b(3)) 
  IF(IZ_max.LT.IZ_min) Return 
     Do I=1,3                         ! constant factorial parameters 
        Do K=1,3 
           J(I+3*K-3)=b(i)-a(k) 
        End do 
     End do 
     J(10)=(j1+j2+j3-3)/2+1 
     Do I=1,3 
        J(I+10)=IZ_min-a(i)               ! initial factorial parameters 
        J(I+13)=b(i)-IZ_min               ! in the sum 
     End do 
     Z=0.d0 
     DO IZ=IZ_min,IZ_max                 ! summation 
        I_max=0                            ! max. factorial 
        Do I=1,16 
           IF (J(i) > I_max) THEN
              I_max=J(i) 
           END IF
        End do 
        Y=1.d0 
        DO I=2,I_max         ! estimation of one term in sum 
           K=0                 ! K - the extent of the integer I in term 
           DO M=1,9 
              IF(J(M) >= I) THEN
                 K=K+1 
              END IF
           End do 
           IF(J(10) >= I) THEN
              K=K-1 
           END IF
           DO M=11,16 
              IF(J(M) >= I) THEN
                 K=K-2 
              END IF
           End do 
           IF(K == 0) Cycle 
              X=DBLE(I)                   ! Y = Y * I ** K/2 
              KK=IABS(K)/2 
              IF(KK > 0 ) THEN 
                 DO M=1,KK 
                    IF(K > 0) THEN
                       Y=Y*X 
                    END IF
                    IF(K < 0) THEN
                       Y=Y/X 
                    END IF
                 END DO 
              END IF 
              IF(mod(K,2) == 1) THEN
                 Y=Y*SQRT(X) 
              END IF
              IF(mod(K,2) == -1) THEN
                 Y=Y/SQRT(X) 
              END IF
        End do 
        IF(mod(IZ,2) == 1) THEN
           Y=-Y 
        END IF
        Z=Z+Y 
        Do I=11,13                  ! new factorial parameters in sum 
           J(I)=J(I)+1 
        End do 
        DO I=14,16 
           J(I)=J(I)-1 
        End do 
     End do                       ! end of summation 
     K=a(1)+a(2)+a(3) 
     IF(mod(k,2) /= 0) THEN
        Z=-Z 
     END IF
     Z_3j=Z 
!***********************************************************************
      END FUNCTION Z_3j 
!***********************************************************************
!***********************************************************************
      FUNCTION Z_6J (j1,j2,j3,j4,j5,j6)
!***********************************************************************
!
!     determination of 6j-symbols without direct using of factorials
!     accoding to formula:
!
!     6j{j1,j2,j3,j4,j5,j6) = {j1,j2,j3}*{j1,j5,j6}*{j4,j2,j3}*{j4,j5,j3}*
!                                SUM(z) {   (-1)^z * (z+1)!   /
!          [ (z-j1-j2-j3)! * (z-j1-j5-j6)! * (z-j4-j2-j3)! *(z-j4-j5-j3)! *
!              (j1+j2+j4+j5-z)! * (j1+j3+j4+j6-z)! * (j2+j3+j5+j6-z)! ]
!
!     where {a,b,c}=sqrt[ (a+b-c)! * (a-b+c)! * (b-a+c)! / (a+b+c+1)! ]
!
!     If we introduce the auxiliary values L(i)
!     (see below the text of program) then
!
!     6j = sqrt{ Pr(j=5,7,i=1,4) (L(j)-L(i))! / Pr(i=1,4) (L(i)+1)! }
!                Sum(z) { (-1)^z * (z+1)! /
!          [ Pr(i=1,4) (z-L(i))!  * Pr(j=5,7) (L(j)-z) ] }
!
!     (the momenta are used in (2J+1)-representation)
!
!----------------------------------------------------------------------

      Implicit none

      Integer,   intent(in)       :: j1,j2,j3,j4,j5,j6
      Integer                     :: I, IZ_min, IZ_max, K, KK, M, IZ, I_max
      Integer                     :: L(7),J(23)
      Real*8                      :: Z_6J, X, R, C

      Z_6J = 0.0

      L(1)=j1+j2+j3-3                    ! auxiliary values
      L(2)=j1+j5+j6-3
      L(3)=j4+j2+j6-3
      L(4)=j4+j5+j3-3
      L(5)=j1+j2+j4+j5-4
      L(6)=j1+j3+j4+j6-4
      L(7)=j2+j3+j5+j6-4
      DO I=1,7
      IF(mod(L(I),2).eq.1) THEN
         Return
      END IF
      L(I)=L(I)/2
      END DO
      IZ_min=MAX0(L(1),L(2),L(3),L(4))   ! limits of the sum
      IZ_max=MIN0(L(5),L(6),L(7))
      IF(IZ_max.LT.IZ_min) THEN
         Return
      END IF
      Do I=1,4
       J(I)=L(I)+1
       Do K=5,7
        M=L(K)-L(I)
        IF(M.LT.0) THEN
           Return                ! check of triangle rule
        END IF
        J(4+3*I+K)=M
       End do
      End do

      Do I=5,8
       J(I)=IZ_min-L(I-4)             ! initial factorial parameters
      End do                          ! in the sum
      Do I=9,11
       J(I)=L(I-4)-IZ_min
      End do
      C=0.0
      DO IZ=IZ_min,IZ_max           ! summation
       I_max=IZ+1
!      this limit for max. factorial follows from symmetry propeties:
!      let's a(i)=L(i) for i=1,4;  b(i)=L(j),j=5,7;
!      then  b(j)-a(i) <= a(k)  <= max[a(k)] = IZ_min;
!      also  a(j) <= b(j), then a(i)-a(j) <= b(i)-a(j) <= IZ_min;
!      and last (a(i)+1) < max(a(i))+1 < IZ_min+1;
       X=1.0
       DO I=2,I_max            ! estimation of one term in sum
        K=2                    ! K - the power of integer I in the term
        DO M=12,23
         IF(J(M).GE.I) THEN
            K=K+1
         END IF
        End do
        DO M=1,4
         IF(J(M).GE.I) THEN
            K=K-1
         END IF
        End do
        DO M=5,11
         IF(J(M).GE.I) THEN
            K=K-2
        END IF
        End do
        IF(K.EQ.0) Cycle
        R=DBLE(I)               ! X = X * I ** K/2
        KK=IABS(K)/2
        IF(KK.GT.0) THEN
         DO M=1,KK
          IF(K.GT.0) X=X*R
          IF(K.LT.0) X=X/R
         END DO
        END IF
        IF(mod(K,2).EQ.+1) THEN
           X=X*SQRT(R)
        END IF
        IF(mod(K,2).EQ.-1) THEN
           X=X/SQRT(R)
        END IF
       End do 
       IF(mod(IZ,2).eq.1) THEN
          X=-X
       END IF
       C=C+X
       Do I=5,8                     ! new factorial parameters in sum
        J(I)=J(I)+1
       End do
       DO I=9,11
        J(I)=J(I)-1
       End do
      End do                        ! end of summation
      Z_6J=C
      END FUNCTION Z_6J
!====================================================================
     FUNCTION Z_6jj(j1,j2,j3,j4,j5,j6)
!====================================================================

      IMPLICIT NONE

      Integer*4, intent(in)             :: j1,j2,j3,j4,j5,j6
      Real*8                            :: Z_6jj
      Real*8, External                   :: Z_6j

      Z_6jj = Z_6j(j1+j1+1,j2+j2+1,j3+j3+1,j4+j4+1,j5+j5+1,j6+j6+1)

      End FUNCTION Z_6JJ
!======================================================================
      FUNCTION Z_9j (j1,j2,j3,j4,j5,j6,j7,j8,j9)
!======================================================================
!
!     determination of 9j-symbols     
!                                    
! {j1 j2 j3}                           {j1 j2 j3} {j4 j5 j6} {j7 j8 j9}        
! {j4 j5 j6} = SUM(j) (-1)^(2j) (2j+1) 
! {j7 j8 j9}                           {j6 j9 j } {j2 j  j8} {j  j1 j4}
!
!
!     (the momenta are used in (2J+1)-representation)
!
!----------------------------------------------------------------------

      Implicit none

      Integer*4, intent(in)                :: j1,j2,j3,j4,j5,j6,j7,j8,j9
      Integer*4                            :: j,i1,i2
      Real*8, External                     :: Z_6j
      Real*8                               :: Z_9j
      
      i1 = max(abs(j1-j9),abs(j2-j6),abs(j4-j8))+1
	  i2 = min(j1+j9,j2+j6,j4+j8)-1

	  Z_9j = 0.d0; if(i1.gt.i2) Return;  if(mod(i2-i1,2).ne.0) Return

      Do j = i1,i2,2
       Z_9j = Z_9j + j * (-1)**(j-1) * &
                     Z_6j(j1,j2,j3,j6,j9,j ) * &
                     Z_6j(j4,j5,j6,j2,j ,j8) * &
                     Z_6j(j7,j8,j9,j ,j1,j4) 
      End do
       
      End FUNCTION Z_9j
END MODULE Angular_Module
