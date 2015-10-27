!************************************************************************************
!deck dvr_input
!***begin prologue     dvr_input
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for dvr basis sets
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       dvr_input
  SUBROUTINE DVR_Input(nphy,nglobal,coord)
  USE dvr_global
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                          :: nphy, nglobal
  CHARACTER(LEN=80)                :: coord
  CHARACTER(LEN=80)                :: chrkey, cprnt
  INTEGER                          :: intkey, lenth, len
  INTEGER                          :: i
  LOGICAL                          :: dollar, logkey, mprnt, atomic
  REAL*8                           :: fpkey
!
  WRITE(iout,1) coord
  LEN=lenth(coord)
!
!     read the input by stopping at the keyword in the input stream.
!
  keywrd='$h0('//coord(1:LEN)//')'
  IF ( dollar(keywrd,card,cpass,inp) ) THEN
!
!      set the print variables for printing.
!
      prloc= 'print='//prnkey
      prloc(12)=chrkey(card,'print',prnkey(12),' ')
      IF(prloc(12) == 'all') THEN
         prn=.true.
      ELSE
         CALL setlog(prn,prloc,card,11)
      END IF
      prloc(12)=chrkey(card,'sector-print',secprn,' ')
      prn(12)=.false.
      IF(prloc(12) == 'sector-details') then
         prn(12)=.true.
      END IF
      cprnt=chrkey(card,'print=dvr_input','on',' ')
      mprnt=.false.
      IF(cprnt == 'on') then
         mprnt = .true.
      END IF
!
!     branch for input depending on whether this is a time 
!     or space axisinate
!
      IF(coord(1:1) /= 't') THEN
!
!       this is a space calculation
!
        WRITE(iout,2)
        mass=fpkey(card,'mass',massau,' ')
        units=chrkey(card,'units','atomic-units',' ')
        parity=chrkey(card,'parity','none',' ')
        angmom=intkey(card,'angular-momentum',0,' ')
        reuse=logkey(card,'reuse-space-data',.false.,' ')
        IF(units == 'atomic-units') THEN
           mass=mass/massau
        END IF
        nodiag=logkey(card,'do-not-diagonalize',.false.,' ')
        typwt=chrkey(card,'weight-type','one',' ')
        atomic=logkey(card,'atomic-calculation',.false.,' ')
        IF(typwt == 'fourier') then
            CALL Read_Fourier
            WRITE(iout,3) edge(1), edge(2)
            WRITE(iout,4) (npt(i),i=1,nreg)
        ELSE IF (typwt == 'legendre') then
            CALL Read_Legendre
            WRITE(iout,4) (npt(i),i=1,nreg)
        ELSE IF (typwt == 'theta') then
            CALL Read_Theta
            WRITE(iout,4) (npt(i),i=1,nreg)
        ELSE IF (typwt == 'hermite') then
            WRITE(iout,5)
            CAll Read_Hermite
        ELSE IF(typwt == 'laguerre') then
            WRITE(iout,6)
            Call Read_laguerre
        ELSE
            CALL Read_Grid_Parameters
        END IF
        write(iout,7) nfix, fix, drop
!
!     Determine nphy
!
        call ptcal(nphy,nglobal,typwt)
        write(iout,8) nphy
!
      ELSE IF (coord(1:1) == 't') THEN
!
!        this is a time calculation
!    
        WRITE(iout,9)
        CALL fparr(card,'region-boundaries',edge,2,' ')
        endpts(1)=edge(1)
        endpts(2)=edge(2)
        n(1)=intkey(card,'polynomial-order',1,' ')
        typwt='one'
        reuse=logkey(card,'reuse-time-data',.false.,' ')
        WRITE(iout,4) (n(i),i=1,2)
        WRITE(iout,10) (edge(i),i=1,2)
!    
!       number of quadrature points is one greater
!       than polynomial order.
    
        npt(1)=n(1)+1
        nphy=n(1)
      ELSE
        write(iout,11)
        call lnkerr('quit:keyword error') 
      END IF
  END IF
1    FORMAT(/,20X,'axis = ',a24)
2    FORMAT(/,20X,'generating dvr representation in space')
3    FORMAT(/,1x,'     Fourier : Region is (',e15.8,',',e15.8,')')
4    FORMAT(/,15x,'polynomial order                       = ',  &
                   (/,15x,5(i4,1x)))                            
5    FORMAT(/,1x,'     Hermite : Region is - infinity to + infinity')
6    FORMAT(/,1x,'     Laguerre: Region is zero to + infinity')
7    Format(/,1x,'number of fixed points = ',i1,/15x,               &
                 'left point fixed       = ',l1,/,15x,              &
                 'right point fixed      = ',l1,/,15x,              &
                 'left point dropped     = ',l1,/15x,               &
                 'right point dropped    = ',l1)
8    FORMAT(/,1x,'     size of the physical basis set = ',i5)
9    FORMAT(/,20X,'generating dvr representation in time')
10   FORMAT(/,15x,'region boundaries                      = ',  &
                   (/,15x,5(e15.8,1x))) 
11   FORMAT(/,1x,'no input found corresponding to keyword = ',a80)

  END SUBROUTINE DVR_Input
!*************************************************************************************
