!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        MODULE normalize_module
                        USE dvrprop_global
!**begin prologue     normalize
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       normalize
!
                        INTERFACE check_norm
             MODULE PROCEDURE check_norm_d,                             &
                              check_norm_z
                    END INTERFACE check_norm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                           CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck check_norm_d.f
!***begin prologue     check_norm_d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            check normalization for three dimension wavefunction
!***
!***references
!***routines called    ddot
!***end prologue       check_norm_d
  SUBROUTINE check_norm_d(v,norm)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                           :: v
  REAL*8                                         :: ddot
  REAL*8                                         :: norm
!
  norm = ddot(n3d,v,1,v,1)
END SUBROUTINE check_norm_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck check_norm_z.f
!***begin prologue     check_norm_z
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            check normalization for three dimension wavefunction
!***
!***references
!***routines called    ddot
!***end prologue       check_norm_z
  SUBROUTINE check_norm_z(v,norm)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)                       :: v
  COMPLEX*16                                     :: cdotc
  REAL*8                                         :: norm
!
  norm = cdotc(n3d,v,1,v,1)
END SUBROUTINE check_norm_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE normalize_module
