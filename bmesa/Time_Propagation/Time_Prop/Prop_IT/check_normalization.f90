!***********************************************************************
                           MODULE check_normalization
                           USE arnoldi_global
!***********************************************************************
                           INTERFACE chk_nrm
                    MODULE PROCEDURE chk_nrm_d, chk_nrm_z
                       END INTERFACE chk_nrm
!***********************************************************************
                           CONTAINS
!***********************************************************************
!***********************************************************************
!***********************************************************************
!deck chk_nrm_d
!**begin prologue     chk_nrm_d
!**date written       010828   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           
!**author             schneider, barry (nsf)
!**source             
!**purpose            check normalization.
!**description        
!**                   
!**references
!**routines called
!**end prologue       chk_nrm_d
  SUBROUTINE chk_nrm_d(v,n)
  IMPLICIT NONE
  INTEGER                     :: n
  REAL*8, DIMENSION(n)        :: v
  REAL*8                      :: ddot
  real*8                      :: nrm
  nrm = ddot(n,v,1,v,1)
  write(iout,1) nrm
1 FORMAT(/,5x,'Normalization Integral = ',e15.8)
END SUBROUTINE chk_nrm_d
!***********************************************************************
!***********************************************************************
!deck chk_nrm_z
!**begin prologue     chk_nrm_z
!**date written       010828   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           
!**author             schneider, barry (nsf)
!**source             
!**purpose            Module to check normalization.
!**description        
!**                   
!**references
!**routines called
!**end prologue       chk_nrm_z
  SUBROUTINE chk_nrm_z(v,n)
  IMPLICIT NONE
  INTEGER                     :: n
  COMPLEX*16, DIMENSION(n)    :: v
  COMPLEX*16                  :: cdotc
  real*8                      :: nrm
  nrm = cdotc(n,v,1,v,1)
  write(iout,1) nrm
1 FORMAT(/,5x,'Normalization Integral = ',e15.8)
END SUBROUTINE chk_nrm_z
!***********************************************************************
!***********************************************************************
END MODULE check_normalization


