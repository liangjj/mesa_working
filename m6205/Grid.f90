!deck Grid.f
!***begin prologue     Grid
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           gnquad, link m6200
!***author             schneider, barry (nsf)
!***source             m6200
!***purpose            Numerical integration quadrature points and weights
!***description        A three dimension quadrature grid is computed using
!***                   either standard product rules in (r,theta,phi) or
!***                   Lebedev quadratures in omega=(theta,phi) and r.
!***references

!***routines called    satshl, shells, voronoi
!***end prologue       Grid
  PROGRAM Grid
  USE Data
  USE Grid_Defined_Types
  IMPLICIT NONE
  INTEGER                                   :: intkey
  INTEGER                                   :: i
  INTEGER                                   :: j
  INTEGER                                   :: iii
  INTEGER                                   :: ibeg
  INTEGER                                   :: iend
  INTEGER                                   :: lenth
  REAL(idp)                                 :: fpkey
  LOGICAL                                   :: logkey
  LOGICAL                                   :: posinp
  LOGICAL                                   :: dollar
  CHARACTER(LEN=3)                          :: itoc

pointer (pz,z(1)), (pz,iz(1))
!

  CALL drum !    Open the input and output files
  CALL iosys ('read character options from rwf',-1,0,0,ops)  ! Read in the options string and the Variables
  Call Set_General_Parameters
  niter=intkey(ops,'grid=number-of-voronoi-iterations',3,' ')
  alpha=fpkey(ops,'grid=exponential-scattering-grid-cutoff', 1.d0,' ')
  prnt=logkey(ops,'grid=print',.false.,' ')
  no_voroni=logkey(ops,'grid=voronoi=off',.false.,' ')
  cutoff=fpkey(ops,'grid=scattering-grid-cutoff',10.d0,' ')
  no_scat=logkey(ops,'grid=no-scattering-center',.false.,' ')  ! Default is scattering center present
  one_grid=logkey(ops,'grid=only-one-grid',.false.,' ')         ! Default is multiple grids
  IF (one_grid) THEN  ! If only one grid, that grid is the scattering grid centered at (0,.0.,0.)
     no_scat=.false.
     no_voroni=.true. ! No Voronoi transformation
  END IF
  no_disk=logkey(ops,'grid=no-disk-output',.false.,' ')
  yukawa_on=logkey(ops,'grid=yukawa=on',.false.,' ')
  yn='yes'
  IF (no_scat==.true.) THEN
    yn='true'
  END IF
  WRITE (iout,1)
  CALL iosys ('read character "grid filename" from rwf',-1,0,0,grid_filename) ! Open the grid file
  CALL iosys ('open grid as new',0,0,0,grid_filename)
  CALL iosys('write character "scattering center" to grid',0,0, 0,yn)  ! Start writing data to the grid file
!            read centers and the parameters for the yukawa potential
  IF ( dollar('$centers',card,cpass,inp) ) THEN
       CALL cardin(card)
  ELSE
       Call lnkerr('error in centers input')
  END IF
  ncent=intkey(card,'no-atomic-centers',1,' ')
  ncplus=ncent
  IF (no_scat==.false.) THEN
      ncplus = ncent + 1
  END IF
  CALL iosys ('write integer "number of atomic centers" to grid',1,ncent,0,' ')
  CALL iosys ('create real "atomic center positions" on grid',3*ncent,0,0,' ')
  CALL iosys ('create real "yukawa exponents" on grid',ncent,0,0,' ')
  CALL iosys ('create real "nuclear charges" on grid',ncent,0,0,' ')
  ALLOCATE(atom(1:ncplus))
  DO i=1,ncent
     ALLOCATE(atom(i)%a(1:3))
     atom(i)%eta=fpkey(card,'exponent-center-'//itoc(i),0.d+00,' ')
     atom(i)%znuc=fpkey(card,'charge-center-'//itoc(i),1.d+00,' ')
     CALL fparr(card,'position-center-'//itoc(i),atom(i)%a, 3,' ')
     WRITE(iout,2) atom(i)%eta, atom(i)%a(1:3), atom(i)%znuc
     CALL iosys ('write real "atomic atomer positions" to grid without rewinding',3,atom(i)%a,0,' ')
     CALL iosys ('write real "yukawa exponents" to grid without rewinding',1,atom(i)%eta,0,' ')
     CALL iosys ('write real "nuclear charges" to grid without rewinding',1,atom(i)%znuc,0,' ')
  END DO
!            scattering center coordinates are read in here. the default is (0.,0.,0.)
  IF (no_scat==.false.) THEN
      CALL fparr(card,'scattering-center-position',atom(ncplus)%a,3,' ')
      WRITE(iout,3) atom(ncplus)%a(1:3)
      CALL iosys ('write real "scattering center position" to grid',3,cent(ncplus)%a,0,' ')
  END IF
!
!                 The following quantities are computed in Satsh in order to
!                 allocate memory.
!
  grd%nrmax=0  ! Largest number of radial points in a single shell needed for any atomic grid
  grd%nthmax=0 ! Largest number of theta points needed for any atomic grid
  grd%nphmax=0 ! Largest number of phi points for needed any atomic grid
  grd%maxgrd=0 ! Largest 3D grid needed for any atomic grid
  grd%maxshl=0 ! Largest number of radial shells needed for any atom
  grd%ntotal=0 ! The total number of words needed to hold all 3D grids
  grd%nwttot=0 ! The total number of words needed to hold all 3D weights
  grd%ltop=0   ! Largest l value for any atom
  grd%mtop=0   ! Largest m value for any atom
  grd%ntrad=0  ! Largest number of radial points needed for any atom
  grd%maxang=0 ! Largest number of angular points for any atom
  ibeg=1
  iend=ncplus
!
!                 Read in the grid information
!
  IF (one_grid) THEN
      ibeg=ncplus
      iend=ncplus
  END IF
  DO  i=ibeg,iend
      IF (i <= ncent) THEN
          cpass='atom'
          WRITE(iout,*)
          WRITE(iout,*) '     reading and computing grid data for atom = ',i
          str='atom-'//itoc(i)
      ELSE
         WRITE(iout,*)
         WRITE(iout,*) '     reading and computing grid data for scattering center'
         cpass='scattering'
         str='scattering-center'
     END IF
!
!    It is important to note that since the defined type grd(i) is passed into
!    the subroutine satshl any variables associates with that grd(i) are available
!    in later calls.  This is an important point when the subroutine shells is called.
!
     call satshl(atom(i))
     atom(i)%ngrid=0    ! total number of points for this atom.  All points
     atom(i)%nwts=0     ! total number of weights for this atom.  If NC radial quadrature its different.
     atom(i)%rad_grd%ntrad=0    ! total number radial points for this atom.
     atom(i)%rad_grd%ntradwt=0  ! total number radial weights for this atom. If NC radial quadrature its different.
     IF (atom(i)%rad_grd%rtyp == 'newton-cotes') THEN
         DO j=1,atom(i)%rad_grd%nshell
            atom(i)%ngrid=atom(i)%ngrid+atom(i)%rad_grd%nr(j)*atom(i)%ang_grd%nang
            grd%nrmax=MAX(grd%nrmax,atom(i)%rad_grd%nr(j))  ! Largest number of radial points for any atomic grid
            atom(i)%rad_grd%nrwts(j)=atom(i)%rad_grd%nr(j)*(atom(i)%rad_grd%nr(j)-1)
            atom(i)%rad_grd%ntrad=atom(i)%rad_grd%ntrad+atom(i)%rad_grd%nr(j)
            atom(i)%rad_grd%ntradwt= atom(i)%rad_grd%ntradwt + atom(i)%rad_grd%nrwts(j)
            atom(i)%nwts=atom(i)%nwts+atom(i)%rad_grd%nrwts(j)*atom(i)%ang_grd%nang
         END DO
     ELSE
         DO j=1,atom(i)%rad_grd%nshell
            atom(i)%ngrid=atom(i)%ngrid+atom(i)%rad_grd%nr(j)*atom(i)%ang_grd%nang
            grd%nrmax=MAX(grd%nrmax,atom(i)%rad_grd%nr(j))
            atom(i)%rad_grd%nrwts(j)=atom(i)%rad_grd%nr(j)
            atom(i)%rad_grd%ntrad=atom(i)%rad_grd%ntrad+atom(i)%rad_grd%nr(j)
            atom(i)%rad_grd%ntradwt= atom(i)%rad_grd%ntradwt + atom(i)%rad_grd%nrwts(j)
            atom(i)%nwts= atom(i)%nwts + atom(i)%rad_grd%nrwts(j)*atom(i)%ang_grd%nang
         END DO
     END IF
     grd%maxgrd=MAX(atom(i)%ngrid,grd%maxgrd)
     grd%ntrad=MAX(atom(i)%rad_grd%ntrad,grd%ntrad)
     grd%ntotal = grd%ntotal + atom(i)%ngrid
     grd%nwttot = grd%nwttot + atom(i)%nwts
     WRITE(iout,4) atom(i)%ngrid
     CALL iosys ('write integer "number of shells '//str//'" to grid',1,atom(i)%rad_grd%nshell,0,' ')
     CALL iosys ('write integer "number radial points per shell '//str//'" to grid',atom(i)%rad_grd%nshell,    &
                  atom(i)%rad_grd%nr,0,' ')
     CALL iosys ('write real "radial edges '//str//'" to grid', atom(i)%rad_grd%nshell+1,atom(i)%rad_grd%r,0,' ')
     CALL iosys('write character "radial quadrature type '//str//'" to grid',0,0,0,atom(i)%rad_grd%rtyp)
     CALL iosys ('write integer "total number of radial points '//str//'" to grid',1,atom(i)%rad_grd%ntrad,0,' ')
     CALL iosys ('write integer "total number of points '//str//'" to grid',1,atom(i)%ngrid,0,' ')
     CALL iosys ('write integer "total number of weights '//str//'" to grid',1,atom(i)%nwts,0,' ')
  END DO
  CALL iosys ('write integer "max number of shells" to grid',1,grd%maxshl,0,' ')
  CALL iosys ('write integer "max l value" to grid',1,grd%ltop,0,' ')
  CALL iosys ('write integer "max m value" to grid',1,grd%mtop,0,' ')
  CALL iosys('write integer "max radial points in shell" to  grid',1,grd%nrmax,0,' ')
  CALL iosys('write integer "max radial points" to grid',1, grd%ntrad,0,' ')
  CALL iosys ('write integer "max theta points" to grid',1, grd%nthetmax,0,' ')
  CALL iosys ('write integer "max phi points" to grid',1, grd%nphmax,0,' ')
  CALL iosys ('write integer "max grid points" to grid',1, grd%maxgrd,0,' ')
  CALL iosys ('write integer "biggest l" to grid',1, grd%ltop,0,' ')
  CALL iosys ('write integer "biggest m" to grid',1, grd%mtop,0,' ')
  IF (nonsep) THEN
      CALL iosys ('write integer "max number lebedev points" to grid',1,grd%maxang,0,' ')
  END IF
  WRITE(iout,5) grd%nrmax, grd%nthmax, grd%nphmax

  ALLOCATE(full_cartesian_grid(1:3,1:ntotal),          &
           full_spherical_grid(1:3,1:ntotal),          &
           full_weights(1:nwttot))
  j1=1
  j2=j1+nrmax*maxshl   !rpt
  j3=j2+MAX(nthmax,maxgrd) !thpt
  j4=j3+nphmax          !phpt
  jadd=j4+nrmax*(nrmax-1)*maxshl   !wtr
  j5=jadd+nrmax     !wtsum
  j6=j5+nthmax      !wtth
  j7=j6+nphmax      !wtph 
  j8=j7+nthmax      !angleb
  j9=j8+nphmax      !work
  j10=j9+nphmax     !wtleb 
  j11=j10+maxgrd
  j12=j11+MAX(nrmax,nthmax,nphmax,ncplus)
  j13=j12+ncplus*ncplus
  j14=j13+3*ntotal
  j15=j14+nwttot
  j16=j15+3*ntotal
  j17=j16+3*maxang
  jhold=j13
  khold=j14
  lhold=j15
!
!      Generate the grids
!
  DO i=ibeg,iend
     IF (i <= ncent) THEN
         cpass='atom'
         chra=itoc(i)
         lstrng=lenth(chra)
         str='atom-'//chra(1:lstrng)
         Write(iout,*)
         Write(iout,*) 'Generating grid for atom = '//str
     ELSE
         str='scattering center'
         cpass='scattering'
         Write(iout,*) 'Generating grid for scattering center'
     END IF
!
!    we are actually storing all of the atomic grids since we
!    need them for the voroni transformation
!
!
  
     ALLOCATE(rpt(1:nrmax,1:maxshl),thpt(1:max(nthmax,maxgrd)),phpt(1:nphmax),                  &
              wtr(1:nrmax*(nrmax-1)*maxshl),wtsum(1:nrmax), wtth(1:nthmax), wtph(1:nphmax),     &
              angleb(1:3*maxang), work(1:max(maxgrd,nrmax*nrmax*maxang)), wtleb(1:maxang),      &
              sthet(1:nthmax), sphi(1:nphmax), cphi(1:nphmax),yuk(1:maxgrd),                    &
              scr(1:max(nrmax,nthmax,nphmax,ncplus)), grid(1:3*ntotal), wt(1:nwttot),           &
              spgrid(1:3*ntotal)

     Call shells(i,rpt,thpt,phpt,wtr,wtsum,wtth,wtph,angleb,work,wtleb,sthet,sphi,cphi,    &
              yuk,scr,grid,wt,spgrid,ns,nrmax,nthet,nphi,ngrid,nwts,ncen,ncplus,prnt,   &
              nleb,nang,nonsep,yukon,nodisk)



  CALL shells(i,z(j1),z(j2),z(j3),z(j4),                &
      z(jadd),z(j5),z(j6),z(j16),z(j17),z(j5),z(j7),    &
      z(j8),z(j9),z(j10),z(j11),z(j13),z(j14),z(j15),   &
      nshell,nrmax,nthet,nphi,ngrid(i),nwts(i),ncent,   &
      ncplus,prnt,nleb,nang,nonsep,yukon,nodisk)
  j13=j13+3*ngrid(i)
  j14=j14+nwts(i)
  j15=j15+3*ngrid(i)
END DO







!        Rather than 
  ALLOCATE(rpt(1:nrmax,1:maxshl),thpt(1:max(nthmax,maxgrd)),phpt(1:nphmax),                  &
           wtr(1:nrmax*(nrmax-1)*maxshl),wtsum(1:nrmax), wtth(1:nthmax), wtph(1:nphmax),     &
           angleb(1:3*maxang), work(1:max(maxgrd,nrmax*nrmax*maxang)), wtleb(1:maxang),      &
           sthet(1:nthmax), sphi(1:nphmax), cphi(1:nphmax),yuk(1:maxgrd),                    &
           scr(1:max(nrmax,nthmax,nphmax,ncplus)), grid(1:3*ntotal), wt(1:nwttot),           &
           spgrid(1:3*ntotal)

  Call shells(i,rpt,thpt,phpt,wtr,wtsum,wtth,wtph,angleb,work,wtleb,sthet,sphi,cphi,    &
              yuk,scr,grid,wt,spgrid,ns,nrmax,nthet,nphi,ngrid,nwts,ncen,ncplus,prnt,   &
              nleb,nang,nonsep,yukon,nodisk)



  CALL shells(i,z(j1),z(j2),z(j3),z(j4),                &
      z(jadd),z(j5),z(j6),z(j16),z(j17),z(j5),z(j7),    &
      z(j8),z(j9),z(j10),z(j11),z(j13),z(j14),z(j15),   &
      nshell,nrmax,nthet,nphi,ngrid(i),nwts(i),ncent,   &
      ncplus,prnt,nleb,nang,nonsep,yukon,nodisk)
  j13=j13+3*ngrid(i)
  j14=j14+nwts(i)
  j15=j15+3*ngrid(i)
END DO
j13=jhold
j14=khold
  IF (.NOT.novoron) THEN
       WRITE(iout,*) 'compute voronoi weights'
       CALL voronoi (eta,a,znuc,z(j13),z(j14),z(j10),z(j1),z(j12),  &
                     z(j11),z(j2),alpha,cutoff,z(j17),ncent,ncplus,  &
                     ngrid,nwts,niter,nr,nrmax,maxshl,noscat, nonsep,prnt,yukon,nodisk)
  ELSE
       WRITE(iout,*) 'do not compute voronoi weights'
END IF
CALL iosys('rewind all on grid read-and-write',0,0,0,' ')
CALL memory(-ngot,pz,idum,'grid',idum)
CALL chainx(0)
1 FORMAT (//,15X,'***** grid generation program *****',//)
2 FORMAT(' yukawa exponent = ',f10.5,' center = (',f10.5,',',f10.5,  &
         ',',f10.5')'/,' charge = ',f10.5)
3 FORMAT(' scattering center ',3F10.5)
4 FORMAT('total number of grid points',1X,i5)
5 FORMAT(/,'maximum radial points in a shell ',i4,  &
        1X,'maximum theta points ',i4,/, 'maximum phi points ',i4)

END PROGRAM Grid
