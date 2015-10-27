!***********************************************************************
                           MODULE Exact_Propagator_Module
                           USE dvrprop_global
                           USE CC_Prop_Module
                           USE initial_state_module
                           USE Iterative_Global
                           USE hamiltonian_module
                           IMPLICIT NONE

!***********************************************************************
!***********************************************************************
                           INTERFACE propagator
             MODULE PROCEDURE propagator_d,                             &
                              propagator_z
                       END INTERFACE propagator
!***********************************************************************
!***********************************************************************
                           Contains
!***********************************************************************
!***********************************************************************
!deck Propagator
!***begin prologue     propagator_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***
!***references
!***routines called
!***end prologue       propagator
!
  SUBROUTINE propagator_d(v_in,v_out)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                     :: v_in, v_out 
!
!
  call ebc(v_out,overlap_d,v_in,n3d,n3d,1)
  call ebtc(v_in,eigenvectors_d,v_out,n3d,n3d,1)
  v_in(:) = EXP(-eigenvalues(:) * deltat/hbar) * v_in(:)
  call ebc(v_out,eigenvectors_d,v_in,n3d,n3d,1)
  title='Exact Value of Exp(-H*t) on Initial Vector'
  call prntfm(title,v_out,n3d,1,n3d,1,iout)
END SUBROUTINE propagator_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck Propagator
!***begin prologue     propagator_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***
!***references
!***routines called
!***end prologue       propagator
!
  SUBROUTINE propagator_z(v_in,v_out)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)      :: v_in, v_out 
!
!
  call cebc(v_out,overlap_z,v_in,n3d,n3d,1)
  call cehbtc(v_in,eigenvectors_z,v_out,n3d,n3d,1)
  v_in(:) = EXP(-eye*eigenvalues(:) * deltat/hbar) * v_in(:)
  call cebc(v_out,eigenvectors_z,v_in,n3d,n3d,1)
  title='Exact Value of Exp(-i*H*t) on Initial Vector'
  call prntcmn(title,v_out,n3d,1,n3d,1,iout,'e')
END SUBROUTINE propagator_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END  MODULE Exact_Propagator_Module

