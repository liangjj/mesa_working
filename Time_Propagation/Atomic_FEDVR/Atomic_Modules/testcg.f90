!
   Program Clebsch_Gordan
!
!  Test of Oleg's Clebsch-Gordan routine
!
!  In this driver, the angular momenta are read in as
!  INTEGERS, with TWICE their actual values.  
!
!  Example: If you want (1/2,1/2;1/2,-1/2|0,0) = 1/sqrt(2)
!
!           the input is   1   1   1   -1  0  0
!
!
      IMPLICIT NONE
      INTEGER                :: J1,M1,J2,M2,J,M 
      INTEGER                :: i 
      INTEGER                :: input=5, output=6 
      INTEGER                :: Num_CG
      REAL*8                 :: CG, CG_Coef 
      OPEN(input,file='CG.inp',status='old')
      OPEN(output,file='CG.out',status='unknown') 
!      print *, 'input  2*j1, 2*m1, 2*j2, 2*m2, 2*j, 2*m'
!      print *, 'first argument < 0 stops the program'
!      print *
      READ(5,*) Num_CG
      DO i=1,Num_CG
         READ(5,*) J1,M1,J2,M2,J,M 
         WRITE (6,1) J1,M1,J2,M2,J,M
         CG = CG_Coef(J1,M1,J2,M2,J,M) 
         WRITE(6,2) CG
      END DO
1     Format(2x,'J1 = ',i4,1x,'M1 = ',i4,1x,'J2 = ',i4,1x,'M2 = ',i4,  &
             1x,'J = ',i4,1x,'M = 'i4)
2     Format(/,5x,'CG_Coef = ',e15.8) 
   END PROGRAM Clebsch_Gordan
!***********************************************************************
!***********************************************************************
!deck CG_Coef
!***begin prologue     CG_Coef
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input, atomic
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            Clebsch-Gordan Coefficients
!***                   
!***description        Clebsch-Gordan Coeficients
!***references
!***routines called    
!***end prologue       CG_Coef
!***********************************************************************
  FUNCTION CG_Coef(L1,M1,L2,M2,L,M)
  IMPLICIT NONE
  INTEGER                                 :: L1, L2, M1, M2, L ,M
  INTEGER, DIMENSION(6)                   :: idum
  REAL*8                                  :: CG_Fac, CG_Coef
  REAL*8, External                        :: Z_3j
  idum(1) = L1 + L1 + 1
  idum(2) = M1 + M1 + 1
  idum(3) = L2 + L2 + 1
  idum(4) = M2 + M2 + 1
  idum(5) = L + L + 1
  idum(6) = M + M + 1
  CG_Fac=(-1)**( ( idum(1) - idum(3) + idum(6)-1)/2 ) * sqrt( DBLE(idum(5)) )
  CG_Coef = CG_Fac * Z_3j( idum(1),idum(2),idum(3),idum(4),idum(5),- idum(6) + 2 )
!
!***********************************************************************
  END FUNCTION CG_Coef
!***********************************************************************
!-------------------------------------------------------------------- 
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
