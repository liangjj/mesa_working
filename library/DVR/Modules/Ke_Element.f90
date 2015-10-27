!deck Ke_Element.f
!***begin prologue     Ke_Element
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            1. calculate piecewise lobatto dvr functions and
!***                      their one-body matrices
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Ke_Element

  SUBROUTINE Ke_Element(pt_0,pt_n,pt,wt,f,df,ddf,ke,eigv_0,eigvec_0,h,eigv, &
                     eigvec,v,kinetic_energy_type,nphy,nglobal)
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                     :: nphy 
  INTEGER                                     :: nglobal 
  CHARACTER(LEN=*)                            :: kinetic_energy_type
  REAL*8, DIMENSION(:), ALLOCATABLE           :: x
  REAL*8, DIMENSION(:), ALLOCATABLE           :: xwt
  REAL*8, DIMENSION(:,:), ALLOCATABLE         :: px
  REAL*8, DIMENSION(:,:), ALLOCATABLE         :: dpx
  REAL*8, DIMENSION(:,:), ALLOCATABLE         :: ddpx
  REAL*8, DIMENSION(:), ALLOCATABLE           :: norm
  INTEGER, DIMENSION(2)                       :: words
  TYPE(regional_mat), DIMENSION(:),                 &
                      ALLOCATABLE        :: reg_mat
  TYPE regional_mat
       REAL*8, DIMENSION(:), ALLOCATABLE      :: qr
       REAL*8, DIMENSION(:), ALLOCATABLE      :: wtr
       REAL*8, DIMENSION(:,:), ALLOCATABLE    :: pr
       REAL*8, DIMENSION(:,:), ALLOCATABLE    :: dpr
       REAL*8, DIMENSION(:,:), ALLOCATABLE    :: ddpr
       REAL*8, DIMENSION(:,:), ALLOCATABLE    :: ovr
       REAL*8, DIMENSION(:,:), ALLOCATABLE    :: tr
       REAL*8, DIMENSION(:), ALLOCATABLE      :: vr
       REAL*8, DIMENSION(:,:), ALLOCATABLE    :: ovr
       REAL*8, DIMENSION(:,:), ALLOCATABLE    :: ham
  END TYPE regional_mat

  REAL*8, DIMENSION(:), ALLOCATABLE  :: glbg
  REAL*8, DIMENSION(:), ALLOCATABLE  :: grid
  REAL*8, DIMENSION(:), ALLOCATABLE  :: mat

  REAL*8                             :: dscale
  INTEGER                            :: i
  INTEGER                            :: start
  INTEGER                            :: end
  INTEGER                            :: err
  INTEGER                            :: hm
  INTEGER                            :: pm
  INTEGER                            :: vm
  INTEGER                            :: km
  character*80                       :: title
!
!
  ALLOCATE(x(nglobal), xwt(nglobal), px(nglobal,nglobal), dpx(nglobal,global),                 &
           ddpx(nglobal,nglobal), norm(nglobal) )
  ALLOCATE( reg_mat(nreg) )
  DO  i=1,nreg
!
!        calculate the un-normalized sector functions and their
!        derivatives.
!
    ALLOCATE( reg_mat(i)%qr(npt(i)), reg_mat(i)%wtr(npt(i)), reg_mat(i)%pr(npt(i),npt(i)),     &
              reg_mat(i)%dpr(npt(i),npt(i)), reg_mat(i)%ddpr(npt(i),npt(i)),                   &
              reg_mat(i)%tr(npt(i),npt(i)), reg_mat(i)%ovr(npt(i),npt(i) )              
    CALL drvply( reg_mat(i)%qr, reg_mat(i)%wtr, reg_mat(i)%pr,                                 &
                 reg_mat(i)%dpr, reg_mat(i)%ddpr, edge(i), typwt, npt(i), i)
!
!        calculate the overlap, bloch and kinetic energy matrix.
!        this is done over un-normalized functions.  renormalization
!        is necessary when joining different sectors.
!
    CALL ovmat( reg_mat(i)%ovr, reg_mat(i)%pr, reg_mat(i)%pr, reg_mat(i)%wtr, npt(i), i)
!
!    Bloch operator now incorporated into kinetic energy and not
!    separately stored.
!
   IF(kinetic_energy_type == 'spheroidal_angular' ) then
             CALL Spheroidal_KE ( reg_mat(i)%tr, reg_mat(i)%pr, reg_mat(i)%dpr,                       &
                          reg_mat(i)%ddpr, reg_mat(i)%qr, reg_mat(i)%wtr, coord, npt(i), i)
    ELSE IF(kinetic_energy_type == 'spheroidal_radial' ) then
             CALL Spheroidal_KE ( reg_mat(i)%tr, reg_mat(i)%pr, reg_mat(i)%dpr,                       &
                          reg_mat(i)%ddpr, reg_mat(i)%qr, reg_mat(i)%wtr, coord, npt(i), i)
    END IF
    CALL vrmat(mat(vr(i)),grid(qr(i)),coord,dscale, &
               npt(i),i)
  END DO
!
!     calculate the global functions
!----------------------------------------------------------------------
!
!                  x = points  xwt=weights px=polynomials
!                  dpx=first derivatives  ddpx=second derivatives
!
!                  these are all full.  no removal of endpoints.
!----------------------------------------------------------------------
!
! In order to use all of the arrays for the subregions, it is necessary
! to pass the array grid and the pointers to routine conbsis.
! For example, glbg(px) is really the array px(nglobal,nglobal) as far
! as the subroutine is concerned.
!
  CALL conbsis(glbg(x),glbg(xwt),glbg(px),glbg(dpx),glbg(ddpx), &
               glbg(norm),grid,qr,wtr,pr,dpr,ddpr,nglobal)
!     calculate the matrix elements and potentials
  hm = 1
  km = hm  + nglobal*nglobal
  vm = km  + nglobal*nglobal
  pm = vm  + nglobal
  need=pm + nglobal*nglobal
  ALLOCATE(ham(need))
  CALL conv(ham(vm),mat,vr,nglobal)
!
!----------------------------------------------------------------------c
!                  hmat=hamiltonian v=potential
!                  first and last functions not removed.
!----------------------------------------------------------------------c
!
! The same situation holds here.  Conham needs the regional submatrices
! and the pointers provide the first word address of these arrays.  The
! subroutine Conham uses these pieces of memory as ordinary arrays of the
! appropriate dimension for the subregion.
!
  CALL conham(ham(hm),ham(vm),ham(km),ham(pm),glbg(norm),mat,tr, &
              p_momr,nglobal)
!     now make all of the physical quantites and diagonalize the
!                          hamiltonian.

  start=1
  END=nglobal
  IF(bcl == 0) THEN
     start=start+1
  END IF
  IF(bcr == 0) THEN
     END=END-1
  END IF
  CALL hamphy(glbg(x),glbg(xwt),glbg(px),glbg(dpx),glbg(ddpx),          &
              ham(hm),ham(vm),ham(km),ham(pm),pt_0,pt_n,pt,wt,f,df,ddf, &
              p_mom,h,v,ke,srf_prm,nglobal,nphy,start)
  IF(.NOT.nodiag) THEN
      write(iout,1)
      CALL diarep(h,eigv,eigvec,nphy)
      CALL rmamp(eigvec,srf_prm,srf,nphy)
  END IF
  IF(.NOT.nodiag) THEN
     write(iout,2)
     CALL diarep(ke,eigv_0,eigvec_0,nphy)
     CALL rmamp(eigvec_0,srf_prm,srf_0,nphy)
  END IF
  DEALLOCATE(ham)
  DEALLOCATE(grid,mat)
  DEALLOCATE(glbg)
1 FORMAT(/,10x,'Diagonalize Full Hamiltonian')
2 FORMAT(/,10x,'Diagonalize Kinetic Energy Hamiltonian')
END SUBROUTINE lobatto
