!deck atomic_input
!***begin prologue     atomic_input
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input, atomic
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for computing dvr basis sets,
!***                   kinetic energy, one-electron and two electron
!***                   radial integrals.
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       atomic_input

  SUBROUTINE atomic_input(nphy,nglobal)
  USE dvr_global
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                          :: nphy, nglobal
  CHARACTER(LEN=80)                :: chrkey
  CHARACTER(LEN=3)                 :: itoc
  INTEGER                          :: intkey, length, len
  INTEGER                          :: i, j, nply, nsubr, nblock, begin
  LOGICAL                          :: dollar, logkey, automte
  REAL*8                           :: fpkey, boundl, boundr, step

!
!     read the input by stopping at the keyword in the input stream.
!
  keywrd='atomic_dvr'
  WRITE(iout,1) keywrd
  keywrd='$'//keywrd
  IF ( dollar(keywrd,card,cpass,inp) ) THEN
!
!      set the print variables for printing.
!
      prloc = 'print='//prnkey
      prloc(12)=chrkey(card,'print',prnkey(12),' ')
      IF(prloc(12) == 'all') THEN
         prn=.true.
      ELSE
         CALL setlog(prn,prloc,card,11)
      END IF
!
!     this is an atomic calculation
!
      drctv=chrkey(card,'type-calculation','all-integrals',' ')
      if(drctv == 'poisson') then
         dentyp=chrkey(card,'type-density','exponential',' ')
      end if
      mass=massau
      units='atomic-units'
      mass=mass/massau
      parity='none'
      charge=fpkey(card,'atomic-charge',1.d0,' ')
      l_orb_max=intkey(card,'maximum-orbital-angular-momentum',0,' ')
      l_max=l_orb_max+l_orb_max
      WRITE(iout,2) charge, l_orb_max, l_max
      nodiag=.true.
      automte=logkey(card,'automate',.false.,' ')
      typwt=chrkey(card,'weight-type','legendre',' ')
      nfix=intkey(card,'number-of-fixed-points',2,' ')
      IF(nfix /= 0) THEN
         fix(1)=logkey(card,'left-fixed-point',.true.,' ')
         fix(2)=logkey(card,'right-fixed-point',.true.,' ')
         drop(1)=logkey(card,'drop-first-function',.false.,' ')
         drop(2)=logkey(card,'drop-last-function',.false.,' ')
      END IF
      bcl=1
      bcr=1
      IF(drop(1)) THEN
         bcl=0
      END IF
      IF(drop(2)) THEN
         bcr=0
      END IF
      if(.not.automte) then
         nreg=intkey(card,'number-of-regions',1,' ')
         CALL fparr(card,'region-boundaries',edge,nreg+1,' ')
         CALL intarr(card,'polynomial-order-per-region',n,nreg,' ')
         npt=n+1
         nrq=npt 
         CALL intarr(card,'number-of-reference-quadrature-'//  &
                          'points-per-region',nrq,nreg,' ')
      ELSE
         nreg=0
         write(iout,*)
         write(iout,*)
         write(iout,*) '                   Automated Selection of', &
                       ' Steps   '
         nblock=intkey(card,'number-of-major-blocks',1,' ')       
         DO i=1,nblock 
            IF ( dollar('$block-'//itoc(i),card, &
                 cpass,inp) ) THEN
                 nply=intkey(card,'default-order',10,' ')
                 nsubr=intkey(card,'number-of-subregions',1,' ')
                 boundl=fpkey(card,'left-boundary',0.d0,' ')
                 boundr=fpkey(card,'right-boundary',0.d0,' ')
                 write(iout,6) i, nply, nsubr, boundl, boundr
                 step = ( boundr - boundl ) / nsubr
                 nreg = nreg + 1
                 begin = nreg
                 edge(nreg)=boundl
                 DO j=1,nsubr
                    nreg = nreg + 1
                    edge(nreg) = edge(nreg-1) + step
                 END DO
                 nreg = nreg - 1
                 n(begin:nreg)=nply
            END IF
         END DO
         npt = n + 1
         nrq = npt 
      END IF
      WRITE(iout,3) nreg
      WRITE(iout,4) (edge(i),i=1,nreg+1)
      WRITE(iout,5) (n(i), i=1,nreg)
!
!     Determine nphy
!
      call ptcal(nphy,nglobal)
      write(iout,6) nphy
      nphy_tri=nphy*(nphy+1)/2
  ELSE
      write(iout,7) keywrd
      call lnkerr('quit')
  END IF
1  FORMAT(/,20X,'keyword = ',a24)
2  FORMAT(/,20X,'generating atomic dvr functions and operators', &
          /,35x,'basic data', &
          /,5x,'charge                           = ',e15.8, & 
          /,5x,'maximum orbital angular momentum = ',i3, &
          /,5x,'maximum total angular momentum   = ',i4 )
3  FORMAT(/,1X,'number of regions = ', i3)
4  FORMAT(/,1X,'space edges = ',/,5(e15.5,1X))
5  FORMAT(/,1X,'polynomial order in each region = ',/,5(i4,1X))
6  FORMAT(/,1X,'basis size = ',i4)
7 FORMAT(/,1x,'no input found corresponding to keyword = ',a80)
END SUBROUTINE atomic_input
