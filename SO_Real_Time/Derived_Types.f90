!
  MODULE Derived_Types
  USE accuracy
  USE Data
  IMPLICIT NONE
!
!***********************************************************************
!***********************************************************************
!____________________________________________________________________________________________!
!____________________________________________________________________________________________!

!
  TYPE CN_LENGTH
     CHARACTER(LEN=80)                                 :: length_title
  END TYPE CN_LENGTH
!
  TYPE CN_VELOCITY
     CHARACTER(LEN=80)                                 :: velocity_title
  END TYPE CN_VELOCITY
!
  TYPE CRANK_NICHOLSON
     TYPE(CN_LENGTH)                                   :: cn_len
     TYPE(CN_VELOCITY)                                 :: cn_vel
     Complex(idp), DIMENSION(:), ALLOCATABLE           :: ham_l
     Complex(idp), DIMENSION(:), ALLOCATABLE           :: ham_u
     Complex(idp), DIMENSION(:), ALLOCATABLE           :: ham_d
     Complex(idp), DIMENSION(:), ALLOCATABLE           :: sol
     Complex(idp), DIMENSION(:), ALLOCATABLE           :: lhs_d  
     Complex(idp), DIMENSION(:), ALLOCATABLE           :: lhs_u
     Complex(idp), DIMENSION(:), ALLOCATABLE           :: lhs_l  
     Complex(idp), DIMENSION(:), ALLOCATABLE           :: rhs_d  
     Complex(idp), DIMENSION(:), ALLOCATABLE           :: rhs_u
     Complex(idp), DIMENSION(:), ALLOCATABLE           :: rhs_l  
     Complex(idp), DIMENSION(:), ALLOCATABLE           :: d  
     Complex(idp), DIMENSION(:), ALLOCATABLE           :: d_l
     Complex(idp), DIMENSION(:), ALLOCATABLE           :: d_u  
     REAL(idp),    DIMENSION(:), ALLOCATABLE           :: pt
     REAL(idp),    DIMENSION(:), ALLOCATABLE           :: wt
!

  END TYPE CRANK_NICHOLSON
!
  TYPE(CRANK_NICHOLSON)                                :: cn
!
  TYPE REAL_H
     REAL(idp),       DIMENSION(:),   ALLOCATABLE      :: eigval
     REAL(idp),       DIMENSION(:,:), ALLOCATABLE      :: eigvec
     REAL(idp),       DIMENSION(:),   ALLOCATABLE      :: ham_l
     REAL(idp),       DIMENSION(:),   ALLOCATABLE      :: ham_u
     REAL(idp),       DIMENSION(:),   ALLOCATABLE      :: ham_d
     REAL(idp),       DIMENSION(:),   ALLOCATABLE      :: sol
     REAL(idp),       DIMENSION(:),   ALLOCATABLE      :: sol_v
     REAL(idp),       DIMENSION(:,:), ALLOCATABLE      :: sol_m
     REAL(idp),       DIMENSION(:,:), ALLOCATABLE      :: H_mat
     REAL(idp),       DIMENSION(:,:), ALLOCATABLE      :: H_eigvec
     REAL(idp),       DIMENSION(:),   ALLOCATABLE      :: H_eigval
     REAL(idp),       DIMENSION(:),   ALLOCATABLE      :: scr
  END TYPE REAL_H
!
  TYPE COMPLEX_H
     REAL(idp),       DIMENSION(:),   ALLOCATABLE      :: eigval
     COMPLEX(idp),    DIMENSION(:,:), ALLOCATABLE      :: eigvec
     COMPLEX(idp),    DIMENSION(:),   ALLOCATABLE      :: ham_l
     COMPLEX(idp),    DIMENSION(:),   ALLOCATABLE      :: ham_u
     COMPLEX(idp),    DIMENSION(:),   ALLOCATABLE      :: ham_d
     COMPLEX(idp),    DIMENSION(:),   ALLOCATABLE      :: sol
     COMPLEX(idp),    DIMENSION(:),   ALLOCATABLE      :: scr
  END TYPE COMPLEX_H

  TYPE COMPLEX_PROP
     COMPLEX(idp),    DIMENSION(:,:), ALLOCATABLE      :: prop
  END TYPE COMPLEX_PROP

  TYPE REAL_PROP
     REAL(idp),       DIMENSION(:,:), ALLOCATABLE      :: prop
  END TYPE REAL_PROP
!
  TYPE MATRICES
     TYPE(COMPLEX_H)                                   :: so_z
     TYPE(REAL_H)                                      :: so_d
     TYPE(COMPLEX_PROP)                                :: pr_z
     TYPE(REAL_PROP)                                   :: pr_d
  END TYPE MATRICES
!
  TYPE ENE_HAM
  END TYPE ENE_HAM
!
  TYPE ENE_EXP
  END TYPE ENE_EXP
!
  TYPE(ENE_HAM)                                        :: hamiltonian
!
  TYPE(ENE_EXP)                                        :: exponential
!
  TYPE IMAG_TIME_SO
     TYPE(MATRICES), DIMENSION(:),   ALLOCATABLE       :: mat       
  END TYPE IMAG_TIME_SO
!
  TYPE REAL_TIME_SO_LEN
     TYPE(MATRICES), DIMENSION(:),   ALLOCATABLE       :: mat 
  END TYPE REAL_TIME_SO_LEN
!
  TYPE REAL_TIME_SO_VEL
     TYPE(MATRICES), DIMENSION(:),   ALLOCATABLE       :: mat 
  END TYPE REAL_TIME_SO_VEL
!
  TYPE SPLIT_OPERATOR
     TYPE(REAL_TIME_SO_LEN)                            :: so_len
     TYPE(REAL_TIME_SO_VEL)                            :: so_vel
     TYPE(IMAG_TIME_SO)                                :: so_imag
     TYPE(MATRICES), DIMENSION(:),   ALLOCATABLE       :: mat 
  END TYPE SPLIT_OPERATOR
!
  TYPE(SPLIT_OPERATOR)                                 :: so
!
!********************************************************************************
  END MODULE Derived_Types
!********************************************************************************
