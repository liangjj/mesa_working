!
MODULE mesa_global
!***begin prologue     mesa_global
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           global mesa variables
!***author             schneider, b. i.(nsf)
!***source             Modules
!***purpose            global variables for mesa
!***description        this routine defines the global variables
!***                   and data needed for the mesa code
!
!***references
 
!***routines called
!***end prologue       mesa_global
  USE io
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
  CHARACTER(LEN=4096)                             :: ops
  CHARACTER(LEN=1600)                             :: card
  CHARACTER (LEN=4)                               :: parity
  CHARACTER(LEN=16), ALLOCATABLE, DIMENSION(:)    :: bflabl
  CHARACTER(LEN=80)                               :: cpass, title
  CHARACTER(LEN=32)                               :: xform
  CHARACTER(LEN=80), DIMENSION(15)                :: prnkey
  CHARACTER(LEN=3), DIMENSION(0:3,0:3,0:3)        :: ctype 
  LOGICAL                                         :: drop
  LOGICAL, DIMENSION(15)                          :: prnloc
  INTEGER                                         :: maxnbf=3000
  INTEGER                                         :: idrop
  INTEGER                                         :: pntbuf
  INTEGER                                         :: nbf, oldnbf, nmotot, nnp
  INTEGER, ALLOCATABLE, DIMENSION(:)              :: list
  INTEGER                                         :: nat, nbasis, nprim, ncont
  INTEGER                                         :: ntypes, nbtype, lenxyz
  INTEGER, ALLOCATABLE, DIMENSION(:)              :: nocart, nobf, maxmom, mintyp
  INTEGER, ALLOCATABLE, DIMENSION(:)              :: minmom, nx, ny, nz, index
  INTEGER, ALLOCATABLE, DIMENSION(:,:)            :: ptprim, noprim, nocont, ptcont, start
  INTEGER                                         :: maxprm, mxcont, maxl, maxblk
  INTEGER                                         :: dlen, dolp
  REAL*8, ALLOCATABLE, DIMENSION(:)               :: ex, cont, zan
  REAL*8, ALLOCATABLE, DIMENSION(:,:)             :: coords, c, s, scr
  DATA prnkey /'print=print=mesa_data_transfer=old-basis',                      &
               'print=print=mesa_data_transfer=new-basis',                      &
               'print=print=mesa_data_transfer=transformation-matrix',          &
               'print=print=mesa_data_transfer=optical-potential',              &
               'print=print=mesa_data_transfer=density-matrix',                 &
               'print=print=mesa_data_transfer=overlap',                        &
               'print=print=mesa_data_transfer=1-electron-direct-hamiltonian',  &
               'print=print=mesa_data_transfer=2-electron-direct-hamiltonian',  &
               'print=print=mesa_data_transfer=direct-hamiltonian',             &  
               'print=print=mesa_data_transfer=kinetic-energy',                 &
               'print=print=mesa_data_transfer=potential-energy',               &
               'mesa_data_transfer=diagonalize',                                &
               'mesa_data_transfer=no-optical-potential',                       &
               'mesa_data_transfer=only-one-electron',                          &
               'mesa_data_transfer=complex-optical-potential'/
!
!
END MODULE mesa_global
