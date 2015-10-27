!
MODULE FEDVR_Shared
!***begin prologue     FEDVR_Shared
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           time, finite element dvr, orthogonal polynomial
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            global shared variables for dvr library
!***description        this routine defines the global variables
!***                   and data needed for the dvrlib
!
!***references

!***routines called    
!***end prologue       FEDVR_Shared
  USE accuracy
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!---------------------------------------------------------------------
!                    shared variables used in lots of routines
!---------------------------------------------------------------------
  CHARACTER (LEN=80), DIMENSION(2)        ::          vtyp
  CHARACTER (LEN=80), DIMENSION(30)       ::          coord
  CHARACTER (LEN=128)                     ::          filbec
  CHARACTER (LEN=8)                       ::          typke
  CHARACTER (LEN=24)                      ::          kinetic_energy_type
  CHARACTER (LEN=24)                      ::          type_calculation
  CHARACTER (LEN=80)                      ::          title
  CHARACTER (LEN=80)                      ::          keyword
  CHARACTER (LEN=80)                      ::          units
  CHARACTER (LEN=80)                      ::          typpot
  CHARACTER (LEN=80)                      ::          typwt
  CHARACTER (LEN=80)                      ::          typarg
  CHARACTER (LEN=80)                      ::          dentyp
  CHARACTER (LEN=80)                      ::          typint
  CHARACTER (LEN=80)                      ::          refwt
  CHARACTER (LEN=80)                      ::          cpass
  CHARACTER (LEN=80)                      ::          drctv
  CHARACTER (LEN=2400)                    ::          card
  CHARACTER (LEN=4)                       ::          parity
  CHARACTER (LEN=16)                      ::          range
  INTEGER                                 ::          physical_points
  INTEGER                                 ::          global_points
  INTEGER, DIMENSION(30)                  ::          nphy
  INTEGER, DIMENSION(30)                  ::          nglobal
  INTEGER, DIMENSION(30)                  ::          row
  INTEGER, DIMENSION(30)                  ::          nonz
  INTEGER, PARAMETER                      ::          maxreg=5000
  INTEGER, DIMENSION(maxreg)              ::          n
  INTEGER, DIMENSION(maxreg)              ::          npt
  INTEGER, DIMENSION(maxreg)              ::          tri_reg
  INTEGER, DIMENSION(maxreg)              ::          nrq
  INTEGER                                 ::          spdim
  INTEGER                                 ::          nfix
  INTEGER                                 ::          bcl
  INTEGER                                 ::          bcr
  INTEGER                                 ::          nreg
  INTEGER                                 ::          l_orb_max 
  INTEGER                                 ::          l_max 
  INTEGER                                 ::          m_max 
  INTEGER                                 ::          l_val 
  INTEGER                                 ::          m_val 
  INTEGER                                 ::          nphy_tri 
  INTEGER                                 ::          angmom
  INTEGER, DIMENSION(:,:), ALLOCATABLE    ::          nfun_reg
  INTEGER, DIMENSION(:), ALLOCATABLE      ::          num_reg
  LOGICAL                                 ::          diag
  LOGICAL                                 ::          nlse
  LOGICAL                                 ::          genpts
  LOGICAL                                 ::          reuse
  LOGICAL                                 ::          nodiag
  LOGICAL                                 ::          reuse_sector_information
  LOGICAL                                 ::          atomic
  LOGICAL                                 ::          proj
  LOGICAL, DIMENSION(2)                   ::          drop
  LOGICAL, DIMENSION(2)                   ::          fix
  LOGICAL                                 ::          skip
  REAL*4,  DIMENSION(20)                  ::          del
  REAL*4,  DIMENSION(20)                  ::          tim
  REAL(idp),  DIMENSION(2)                ::          endpts
  REAL(idp)                               ::          mass
  REAL(idp),  DIMENSION(maxreg+1)         ::          edge
  REAL(idp)                               ::          box
  REAL(idp)                               ::          deltax
  REAL(idp)                               ::          dscale
  REAL(idp)                               ::          Z_a
  REAL(idp)                               ::          Z_b
  REAL(idp)                               ::          R_ab
  REAL(idp)                               ::          pre_factor
  REAL(idp)                               ::          factor
  REAL(idp)                               ::          R_max
  LOGICAL                                 ::          toau
  LOGICAL                                 ::          useau
  LOGICAL                                 ::          prnt
!
!
!---------------------------------------------------------------------
!                    Keywords for Print Options
!---------------------------------------------------------------------
!
  CHARACTER (LEN=80), DIMENSION(12) :: prnkey
  CHARACTER (LEN=80)                :: secprn
  DATA prnkey /'sector_points', 'sector_factors',                         &
               'sector_polynomials', 'sector_matrices',                   &
               'global_points', 'global_polynomials',                     &
               'potential','global_matrices',                             &
               'hamiltonian', 'eigenvalues','eigenvectors', 'all' /
  DATA secprn / 'none' /
!
!
 CHARACTER (LEN=80),DIMENSION(12)         :: prloc
 LOGICAL, DIMENSION(12)                   :: prn
!
!
END MODULE FEDVR_Shared
