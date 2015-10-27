!deck lobatto.f
!***begin prologue     lobatto
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
!***end prologue       lobatto

  SUBROUTINE lobatto(pt_0,pt,wt,f,df,ddf,ke,p_mom,eigv_0,eigvec_0,h,eigv, &
                     eigvec,v,srf_prm,srf_0,srf,coord,nphy,nglobal)
  USE dvr_global
  IMPLICIT NONE
  INTEGER                            :: nphy, nglobal 
  REAL*8                             :: pt_0
  REAL*8, DIMENSION(nphy)            :: pt
  REAL*8, DIMENSION(nphy)            :: wt, eigv_0, eigv, v
  REAL*8, DIMENSION(nphy,nphy)       :: f, df, ddf
  REAL*8, DIMENSION(nphy,nphy)       :: ke, p_mom, eigvec_0
  REAL*8, DIMENSION(nphy,nphy)       :: h, eigvec
  REAL*8, DIMENSION(2)               :: srf_prm
  REAL*8, DIMENSION(nphy,2)          :: srf_0, srf
  CHARACTER(LEN=*)                   :: coord
  INTEGER                            :: need, wrd
  INTEGER                            :: x, xwt, px, dpx, ddpx, norm
  INTEGER, DIMENSION(2)              :: words
  INTEGER, DIMENSION(nreg)           :: qr, wtr, pr, dpr, ddpr      
  INTEGER, DIMENSION(nreg)           :: ovr, tr, p_momr, vr               
  REAL*8, DIMENSION(:), ALLOCATABLE  :: glbg, grid, mat, ham
  REAL*8, DIMENSION(:), ALLOCATABLE  :: work
  REAL*8                             :: dscale
  INTEGER                            :: i, start, end, err
  INTEGER                            :: hm, pm, vm, km
  character*80                       :: title
!
INTERFACE 
  Subroutine MEMREG(qr,wtr,pr,dpr,ddpr,ovr,tr,p_momr,vr,words)
  USE dvr_global, ONLY  : nreg, npt  
  IMPLICIT NONE
  INTEGER, DIMENSION(nreg) :: qr, wtr, pr, dpr, ddpr 
  INTEGER, DIMENSION(nreg) :: ovr, tr, p_momr, vr
  INTEGER, dimension(2)    :: words
  END Subroutine MEMREG
END INTERFACE

!
!     Get some memory for arrays used to hold the dvr points, weights,
!     functions, first and second derivatives on the global grid.
!     The integer arrays x, xwt, px, dpx, ddpx, norm, qr, wtr, pr, dpr, ddpr
!     ovr, tr, and p_momr are introduced to point to each array in allocated list.  
!     I have avoided using defined types and pointers to pointers to be able 
!     to allocate a single, large workspace and to partition it out as I need it.
!
  x=1
  xwt=x+nglobal
  px=xwt+nglobal
  dpx=px+nglobal*nglobal
  ddpx=dpx+nglobal*nglobal
  norm=ddpx+nglobal*nglobal
  wrd=norm+nglobal-1
  need=wrd+1
  ALLOCATE(glbg(need),stat=err)
  if(err /= 0 )then
     call lnkerr('allocation error')
  end if
!  CALL rzero(glbg,wrd)
  glbg=0.d0
!
! Memreg returns the first word address of the arrays for each of the regions.
! In subsequent call we pass the allocated array grid using these pointers
! so that they may be treated as ordinary arrays for each region, dimensioned
! as one would expect.
!
  call memreg(qr,wtr,pr,dpr,ddpr,ovr,tr,p_momr,vr,words)
  ALLOCATE(grid(words(1)),mat(words(2)),stat=err)
  if(err /= 0 )then
     call lnkerr('allocation error')
  end if
  DO  i=1,nreg
!
!        calculate the un-normalized sector functions and their
!        derivatives.
    CALL drvply(grid(qr(i)),grid(wtr(i)),grid(pr(i)), &
                grid(dpr(i)),grid(ddpr(i)),edge(i),typwt, &
                parity,angmom,npt(i),nrq(i))
!
!        calculate the overlap, bloch and kinetic energy matrix.
!        this is done over un-normalized functions.  renormalization
!        is necessary when joining different sectors.
!
    CALL ovmat(mat(ovr(i)),grid(pr(i)),grid(pr(i)), &
               grid(wtr(i)),npt(i),i)
!
!    Bloch operator now incorporated into kinetic energy and not
!    separately stored.
!
    IF (typwt == 'legendre') then
        CALL ke_legendre(mat(tr(i)),grid(pr(i)),grid(dpr(i)),     &
                         grid(ddpr(i)),grid(qr(i)),grid(wtr(i)),  &
                         coord,npt(i),i)
    ELSEIF (typwt == 'hermite') then
        CALL ke_hermite(mat(tr(i)),grid(pr(i)),grid(pr(i)), &       
                        grid(dpr(i)),grid(ddpr(i)),grid(qr(i)), &
                        grid(wtr(i)),npt(i),i)
        CALL p_hermite (mat(p_momr(i)),grid(pr(i)),grid(pr(i)), &
                        grid(dpr(i)),grid(qr(i)),grid(wtr(i)),npt(i),i)
    ELSE IF(typwt == 'laguerre' ) then
        CALL ke_laguerre(mat(tr(i)),grid(pr(i)),grid(pr(i)), &       
                         grid(dpr(i)),grid(ddpr(i)),          &
                         grid(wtr(i)),npt(i),i)
        CALL p_laguerre (mat(p_momr(i)),grid(pr(i)),grid(pr(i)), &
                         grid(dpr(i)),grid(wtr(i)),npt(i),i)
    ELSE
        CALL kemat(mat(tr(i)),grid(pr(i)),grid(dpr(i)), &
                   grid(ddpr(i)),grid(qr(i)),grid(wtr(i)),coord,npt(i),i)
        CALL pmat (mat(p_momr(i)),grid(pr(i)),grid(dpr(i)), &
                   grid(qr(i)),grid(wtr(i)),coord,npt(i),i)
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
  CALL hamphy(glbg(x),glbg(xwt),glbg(px),glbg(dpx),glbg(ddpx),  &
              ham(hm),ham(vm),ham(km),ham(pm),pt_0,pt,wt,f,df,ddf, &
              p_mom,h,v,ke,srf_prm,nglobal,nphy,start)
  IF(.NOT.nodiag) THEN
      ALLOCATE(work(1:5*nphy))
      write(output,1)
      CALL diarep(h,eigv,eigvec,work,nphy)
      CALL rmamp(eigvec,srf_prm,srf,nphy)
  END IF
  IF(.NOT.nodiag) THEN
     write(output,2)
     CALL diarep(ke,eigv_0,eigvec_0,work,nphy)
     CALL rmamp(eigvec_0,srf_prm,srf_0,nphy)
  END IF
  IF(.NOT.nodiag) THEN
      DEALLOCATE(work)
  END IF
  DEALLOCATE(ham)
  DEALLOCATE(grid,mat)
  DEALLOCATE(glbg)
1 FORMAT(/,10x,'Diagonalize Full Hamiltonian')
2 FORMAT(/,10x,'Diagonalize Kinetic Energy Hamiltonian')
END SUBROUTINE lobatto
