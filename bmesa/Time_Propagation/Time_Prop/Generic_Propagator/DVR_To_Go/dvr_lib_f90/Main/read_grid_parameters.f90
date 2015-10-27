!************************************************************************************
!deck read_grid_parameters
!***begin prologue     read_grid_parameters
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for dvr basis sets
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       read_grid_parameters
  SUBROUTINE Read_Grid_Parameters
  USE dvr_global
  USE dvr_prnt
  IMPLICIT NONE
  CHARACTER(LEN=3)                 :: itoc
  INTEGER                          :: intkey
  INTEGER                          :: i, j, nply, nsubr, nblock, begin, ntest
  INTEGER                          :: n_sub
  LOGICAL                          :: dollar, logkey, automte
  REAL*8                           :: fpkey, step, boundl, boundr
  automte=logkey(card,'automate',.false.,' ')
  nfix=intkey(card,'number-of-fixed-points',2,' ')
  fix(1)=.false.
  fix(2)=.false.
  drop(1)=.false.
  drop(2)=.false.
  IF(nfix /= 0) THEN
     fix(1)=logkey(card,'left-fixed-point',.false.,' ')
     fix(2)=logkey(card,'right-fixed-point',.false.,' ')
     drop(1)=logkey(card,'drop-left-function',.false.,' ')
     drop(2)=logkey(card,'drop-right-function',.false.,' ')
  END IF
  bcl=1
  bcr=1
  IF(drop(1)) THEN
     bcl=0
  END IF
  IF(drop(2)) THEN
     bcr=0
  END IF
  IF(.not.automte) THEN
     nreg=intkey(card,'number-of-regions',1,' ')
     CALL fparr(card,'region-boundaries',edge,nreg+1,' ')
     CALL intarr(card,'polynomial-order-per-region',n,nreg,' ')
     npt=n+1
     nrq=npt 
     CALL intarr(card,'number-of-reference-quadrature-'//  &
                      'points-per-region',nrq,nreg,' ')
  ELSE
     write(iout,*)
     write(iout,*)
     write(iout,*) '                   Automated Selection of', &
                   ' Steps   '
     nreg=0
     nblock=intkey(card,'number-of-major-blocks',1,' ')       
     DO i=1,nblock 
        IF ( dollar('$block-'//itoc(i),card, &
             cpass,inp) ) THEN
             skip=logkey(card,'skip',.false.,' ')
             IF(skip) THEN
                write(iout,2) i
             ELSE 
                nply=intkey(card,'default-order',10,' ')
                nsubr=intkey(card,'number-of-subregions',1,' ')
                boundl=fpkey(card,'left-boundary',0.d0,' ')
                boundr=fpkey(card,'right-boundary',0.d0,' ')
                n_sub=intkey(card,'number-of-reference-quadrature-'//  &
                                  'points-per-sub-region',2*nply,' ')
                write(iout,1) i, nply, nsubr, boundl, boundr
                step = ( boundr - boundl ) / nsubr
                nreg = nreg + 1
                begin = nreg
                edge(nreg)=boundl
                DO j=1,nsubr
                   nrq(nreg)=n_sub
                   nreg = nreg + 1
                   edge(nreg) = edge(nreg-1) + step
                END DO
                nreg = nreg - 1
                n(begin:nreg)=nply
             END IF
        END IF
     END DO
     npt = n + 1
!     nrq = npt 
  END IF
1    FORMAT(/,1x, 'block = ',i3, &
            /,15x,'polynomial order     = ',i4, &
            /,15x,'number of subregions = ',i4, & 
            /,15x,'left hand boundary   = ',e15.8, &
            /,15x,'right hand boundary  = ',e15.8)
2   FORMAT(/,1x,'skipping input block = ',i4)
  END SUBROUTINE Read_Grid_Parameters
