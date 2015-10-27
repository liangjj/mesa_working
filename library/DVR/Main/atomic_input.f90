!deck spheroidal_input
!***begin prologue     spheroidal_input
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input, spheroidal
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for computing dvr basis sets,
!***                   kinetic energy, one-electron and two electron
!***                   radial integrals.
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       spheroidal_input

  SUBROUTINE spheroidal_input(nphy,nglobal)
  USE dvr_global
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                          :: nphy, nglobal
  CHARACTER(LEN=80)                :: chrkey
  CHARACTER(LEN=3)                 :: itoc
  INTEGER                          :: intkey, lenth, len
  INTEGER                          :: i, j, nply, nsubr, nblock, begin
  LOGICAL                          :: dollar, logkey
  LOGICAL, DIMENSION(2)            :: automte
  REAL*8                           :: fpkey, boundl, boundr, step

!
!     read the input by stopping at the keyword in the input stream.
!
  keywrd='spheroidal_dvr'
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
!     this is an spheroidal calculation
!
      drctv=chrkey(card,'type_calculation','all_integrals',' ')
      if(drctv == 'poisson') then
         dentyp=chrkey(card,'type_density','exponential',' ')
      end if
      mass=massau
      units='atomic_units'
      mass=mass/massau
      parity='none'
      z_a=fpkey(card,'charge_nucleus_a',1.d0,' ')
      z_b=fpkey(card,'charge_nucleus_b',1.d0,' ')
      R_ab=fpkey(card,'internuclear_distance',1.d0,' ')
      l_max=intkey(card,'maximum_l_value',0,' ')
      m_max=intkey(card,'maximum_m_value',l_max,' ')
      WRITE(iout,2) z_a, z_b, R_ab, l_max, m_max
      nodiag=.true.
      typwt=chrkey(card,'weight_type','legendre',' ')
  END IF
  keywrd='$eta_variable'
  IF ( dollar(keywrd,card,cpass,inp) ) THEN
      automte(1)=logkey(card,'automate',.false.,' ')
      nfix(1)=intkey(card,'number_of_fixed_points',1,' ')
      IF(nfix(1) /= 0) THEN
         fix(1,1)=logkey(card,'left_fixed_point',.false.,' ')
         fix(1,2)=logkey(card,'right_fixed_point',.false.,' ')
         drop(1,1)=logkey(card,'drop_first_function',.false.,' ')
         drop(1,2)=logkey(card,'drop_last_function',.false.,' ')
         bcl(1)=1
         bcr(1)=1
         IF(drop(1,1)) THEN
            bcl(1)=0
         END IF
         IF(drop(1,2)) THEN
            bcr(1)=0
         END IF
      END IF
  END IF
  keywrd='$xi_variable'
  IF ( dollar(keywrd,card,cpass,inp) ) THEN
      automte(2)=logkey(card,'automate',.false.,' ')
      nfix(2)=intkey(card,'number_of_fixed_points',1,' ')
      IF(nfix(2) /= 0) THEN
         fix(2,1)=logkey(card,'left_fixed_point',.false.,' ')
         fix(2,2)=logkey(card,'right_fixed_point',.false.,' ')
         drop(2,1)=logkey(card,'drop_first_function',.false.,' ')
         drop(2,2)=logkey(card,'drop_last_function',.false.,' ')
         bcl(2)=1
         bcr(2)=1
         IF(drop(2,1)) THEN
            bcl(2)=0
         END IF
         IF(drop(2,2)) THEN
            bcr(2)=0
         END IF
      END IF
  END IF
      if(.not.automte(1)) then
         nreg=intkey(card,'number_of_regions',1,' ')
         CALL fparr(card,'region_boundaries',edge,nreg+1,' ')
         CALL intarr(card,'polynomial_order_per_region',n,nreg,' ')
         npt=n+1
         nrq=npt 
         CALL intarr(card,'number_of_reference_quadrature_'//  &
                          'points_per_region',nrq,nreg,' ')
      ELSE
         nreg=0
         write(iout,*)
         write(iout,*)
         write(iout,*) '                   Automated Selection of', &
                       ' Steps   '
         nblock=intkey(card,'number_of_major_blocks',1,' ')       
         DO i=1,nblock 
            IF ( dollar('$block_'//itoc(i),card, &
                 cpass,inp) ) THEN
                 nply=intkey(card,'default_order',10,' ')
                 nsubr=intkey(card,'number_of_subregions',1,' ')
                 boundl=fpkey(card,'left_boundary',0.d0,' ')
                 boundr=fpkey(card,'right_boundary',0.d0,' ')
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
1  FORMAT(/,20X,'keyword = ',a24
2  FORMAT(/,20X,'generating atomic dvr functions and operators', &
          /,35x,'basic data', &
          /,5x,'charge on nucleus a              = ',e15.8, & 
          /,5x,'charge on nucleus b              = ',e15.8, & 
          /,5x,'internuclear distance            = ',e15.8, & 
          /,5x,'maximum l value                  = ',i3, &
          /,5x,'maximum m value                  = ',i4 )
3  FORMAT(/,1X,'number of regions = ', i3)
4  FORMAT(/,1X,'space edges = ',/,5(e15.5,1X))
5  FORMAT(/,1X,'polynomial order in each region = ',/,5(i4,1X))
6  FORMAT(/,1X,'basis size = ',i4)
7 FORMAT(/,1x,'no input found corresponding to keyword = ',a80)
END SUBROUTINE spheroidal_input
