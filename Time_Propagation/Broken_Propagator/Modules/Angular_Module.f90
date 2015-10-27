!***********************************************************************
                           MODULE Angular_Module
                           IMPLICIT NONE
!***********************************************************************
!!***********************************************************************
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
END MODULE Angular_Module
