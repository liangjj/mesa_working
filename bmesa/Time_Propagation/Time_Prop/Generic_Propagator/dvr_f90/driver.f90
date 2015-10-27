!deck driver.f
!***begin prologue     driver
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
!***end prologue       driver
  PROGRAM driver
  USE dvr_global
  IMPLICIT NONE
  INTEGER                :: nphy, nglobal, dim, row, intkey, i, err, k
  CHARACTER (LEN=2)      :: itoc
  CHARACTER (LEN=80)     :: chrkey, coord, title
  CHARACTER (LEN=8)      :: typke
  CHARACTER*4096 ops
  LOGICAL                :: dollar
  REAL*8, DIMENSION(4)                   :: pt_0, pt_n(4)
  TYPE dvr_grid
    REAL*8, DIMENSION(:),   POINTER      :: pt, wt, eigv_0, eigv, v
    REAL*8, DIMENSION(:,:), POINTER      :: f, df, ddf, ke, eigvec_0
    REAL*8, DIMENSION(:,:), POINTER      :: eigvec, h, p_mom
    REAL*8, DIMENSION(:),   POINTER      :: srf_prm
    REAL*8, DIMENSION(:,:), POINTER      :: srf_0, srf
  END TYPE dvr_grid
  TYPE (dvr_grid), DIMENSION(:), ALLOCATABLE &
                                         :: grid
  INTEGER                                :: j
!  CALL Drum
  open(inp,file='dvr.inp',status='old')
  open(iout,file='dvr.out',status='unknown')
!  call iosys ('read character options from rwf',-1,0,0,ops) 
  IF( dollar ('$dvr_input',card,cpass,inp ) ) then
      dim=intkey(card,'number-of-space-variables',1,' ')
      typke=chrkey(card,'kinetic-energy-type','dvr',' ')
      do i=1,dim
         coord=chrkey(card,'space-variable-'//itoc(i),'x',' ')      
         write(iout,1) coord
         ALLOCATE(grid(i))
         if(typke == 'dvr'.OR.typke == 'packed') THEN
            call DVR_Input(nphy,nglobal,coord)
            write(iout,2) nphy, nglobal
            IF(.not.nodiag) then
                ALLOCATE(grid(i)%eigv_0(nphy),grid(i)%eigv(nphy), &
                         grid(i)%eigvec_0(nphy,nphy),             &
                         grid(i)%eigvec(nphy,nphy),stat=err)
                if (err /= 0 ) then
                    call lnkerr('allocation error')
                end if
            END IF
            IF(typwt /= 'fourier') THEN
                ALLOCATE(grid(i)%pt(nphy),          &
                         grid(i)%wt(nphy),          &
                         grid(i)%f(nphy,nphy),      &
                         grid(i)%df(nphy,nphy),     &
                         grid(i)%ddf(nphy,nphy),    &
                         grid(i)%ke(nphy,nphy),     &
                         grid(i)%p_mom(nphy,nphy),  &
                         grid(i)%h(nphy,nphy),      &
                         grid(i)%v(nphy),           &
                         grid(i)%srf_prm(2),        &
                         grid(i)%srf_0(nphy,2),     &
                         grid(i)%srf(nphy,2),       &
                         stat=err)
                if (err /= 0 ) then
                    call lnkerr('allocation error')
                end if
                call dvr_basis(pt_0(i),pt_n(i),grid(i)%pt,       &
                               grid(i)%wt,grid(i)%f,grid(i)%df,  &
                               grid(i)%ddf,grid(i)%ke,           &
                               grid(i)%p_mom,grid(i)%eigv_0,     &
                               grid(i)%eigvec_0,grid(i)%h,       &
                               grid(i)%eigv,grid(i)%eigvec,      &
                               grid(i)%v,grid(i)%srf_prm,        &
                               grid(i)%srf_0,grid(i)%srf,coord,  &
                               nphy,nglobal)
                DEALLOCATE(grid(i)%pt,              &
                           grid(i)%wt,              &
                           grid(i)%f,               &
                           grid(i)%df,              &
                           grid(i)%ddf,             &
                           grid(i)%ke,              &
                           grid(i)%p_mom,           &
                           grid(i)%h,               &
                           grid(i)%v,               &
                           grid(i)%srf_prm,         &
                           grid(i)%srf_0,           &
                           grid(i)%srf)
            ELSE
                ALLOCATE(grid(i)%pt(nphy),          &
                         grid(i)%wt(nphy),          &
                         grid(i)%f(nphy,nphy),      &
                         grid(i)%ke(nphy,nphy),     &
                         grid(i)%h(nphy,nphy),      &
                         grid(i)%v(nphy),           &
                         stat=err)
                if (err /= 0 ) then
                    call lnkerr('allocation error')
                end if
                call fourier(grid(i)%pt,grid(i)%wt,  &
                             grid(i)%f,edge,nphy)
                call fourier_basis(grid(i)%pt,       &
                                   grid(i)%ke,       &
                                   grid(i)%h,        &
                                   grid(i)%v,        &
                                   grid(i)%eigv_0,   &
                                   grid(i)%eigvec_0, &
                                   grid(i)%eigv,     &
                                   grid(i)%eigvec,   &
                                   coord,nphy)
                DEALLOCATE(grid(i)%pt,              &
                           grid(i)%wt,              &
                           grid(i)%f,               &
                           grid(i)%ke,              &
                           grid(i)%h,               &
                           grid(i)%v)
            END IF
            IF(.not.nodiag) then
               DEALLOCATE(grid(i)%eigv_0,           &
                          grid(i)%eigv,             &
                          grid(i)%eigvec_0,         &
                          grid(i)%eigvec)
            END IF
         ELSE
               CALL fd_input(nphy,nglobal,row,coord)
               ALLOCATE(grid(i)%pt(nphy), &
                        grid(i)%wt(nphy), &
                        grid(i)%ke(row,nphy), &
                        grid(i)%h(row,nphy), &
                        grid(i)%v(nphy))
               IF(.not.nodiag) then
                  ALLOCATE(grid(i)%eigv_0(nphy), &
                           grid(i)%eigvec_0(nphy,nphy), &
                           grid(i)%eigv(nphy), &
                           grid(i)%eigvec(nphy,nphy))
               endif
               CALL fd_basis(pt_0(i),grid(i)%pt,grid(i)%wt, &
                             grid(i)%ke,grid(i)%eigv_0,grid(i)%eigvec_0, &
                             grid(i)%h,grid(i)%eigv,grid(i)%eigvec, &
                             grid(i)%v,nphy,nglobal,row,coord)
               DEALLOCATE(grid(i)%pt, &
                          grid(i)%wt, &
                          grid(i)%ke, &
                          grid(i)%h, &
                          grid(i)%v)
               IF(.not.nodiag) then
                  DEALLOCATE(grid(i)%eigv_0, &
                             grid(i)%eigvec_0, &
                             grid(i)%eigv, &
                             grid(i)%eigvec)
               endif
         END IF
      END DO
  ELSE
      write(iout,3)
      stop
  END IF
!  CALL Chainx(0)
  stop
1    FORMAT(/,20X,'One-Body DVR Basis Code. coord = ',a24)
2    FORMAT(/,5x,'Nphy = ',i5,2x,'Nglobal = ',i5)
3    FORMAT(/,10x,'Error In Input File')
END PROGRAM driver
