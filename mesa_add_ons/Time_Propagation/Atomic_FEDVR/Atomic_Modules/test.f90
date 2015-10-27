!***********************************************************************
                           MODULE Atomic_Module
                           USE dvr_shared
                           USE dvr_global
                           USE dvrprop_global
                           IMPLICIT NONE
             INTEGER             :: n_tri 
             INTEGER             :: l_orb_max 
             INTEGER             :: l_max 
             CHARACTER(LEN=80)   :: ints, dentyp
             REAL*8              :: atomic_charge
             CHARACTER(LEN=80)   :: chrkey
             CHARACTER(LEN=3)    :: itoc
             LOGICAL             :: dollar, logkey
             INTEGER             :: intkey
             REAL*8              :: fpkey
             REAL*8, DIMENSION(:,:),                             &
                     ALLOCATABLE  :: h_one
             REAL*8, DIMENSION(:,:),                             &
                     ALLOCATABLE  :: v_two
             TYPE Orbitals
               INTEGER, DIMENSION(:), POINTER    :: M_Val
             END TYPE Orbitals
             TYPE (Orbitals), DIMENSION, ALLOCATABLE             &
                                                 :: Orb
             LOGICAL, DIMENSION(:),                              &
                      ALLOCATABLE :: M_Val
END Atomic_Module
