!deck m6204
!***begin prologue     m6204
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           gnquad, link m6200
!***author             schneider, barry (nsf)
!***source             m6204
!***purpose            Numerical integration quadrature points and weights
!***description        A three dimension quadrature grid is computed using
!***                   either standard product rules in (r,theta,phi) or
!***                   Lebedev quadratures in omega=(theta,phi) and r.
!***references

!***routines called    satshl, shells, voronoi
!***end prologue       m6204
  PROGRAM m6204
  USE Data
  USE Grid_Defined_Types
  USE Shell_Info
  IMPLICIT NONE
!  TYPE(ATOMS), DIMENSION(:), ALLOCATABLE    :: atom
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
!
  CALL drum !    Open the input and output files
  CALL iosys ('read character options from rwf',-1,0,0,ops)  ! Read in the options string and the Variables
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
      CALL iosys ('write real "scattering center position" to grid',3,atom(ncplus)%a,0,' ')
  END IF
!
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

!    in later calls.  This is an important point when the subroutine shells is called.
!
     call satshl(atom(i))
  END DO
  CALL iosys('rewind all on grid read-and-write',0,0,0,' ')
  CALL chainx(0)
1 FORMAT (//,15X,'***** grid generation program *****',//)
2 FORMAT(' yukawa exponent = ',f10.5,' center = (',f10.5,',',f10.5,  &
         ',',f10.5')'/,' charge = ',f10.5)
3 FORMAT(' scattering center ',3F10.5)
4 FORMAT('total number of grid points',1X,i5)
5 FORMAT(/,'maximum radial points in a shell ',i4,  &
        1X,'maximum theta points ',i4,/, 'maximum phi points ',i4)

END PROGRAM m6204
