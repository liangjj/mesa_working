!deck dvr_function.f
!***begin prologue     dvr_function
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr_function functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Compute the DVR functions and 
!***                   first and second derivatives
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       dvr_function

  SUBROUTINE dvr_function(f,df,ddf,p,dp,ddp,wtfn,dwtfn,ddwtfn,pt,wt,nphy)
  USE dvr_global
  IMPLICIT NONE
  INTEGER                            :: nphy 
  REAL*8, DIMENSION(nphy)            :: pt, wt
  REAL*8, DIMENSION(nphy,nphy)       :: f, df, ddf
  REAL*8, DIMENSION(nphy,nphy)       :: p, dp, ddp
  REAL*8, DIMENSION(nphy)            :: wtfn, dwtfn, ddwtfn
  character*80                       :: title
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
      write(iout,1)
      CALL diarep(h,eigv,eigvec,work,nphy)
      CALL rmamp(eigvec,srf_prm,srf,nphy)
  END IF
  IF(.NOT.nodiag) THEN
     write(iout,2)
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
END SUBROUTINE dvr_function
