!dvr_input.f90,v 1.2 2004/11/22 18:57:15 bis Exp bis $
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

  SUBROUTINE dvr_input(nphy,nglobal,coord)
  USE dvr_global
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                          :: nphy, nglobal
  CHARACTER(LEN=*)                 :: coord
  CHARACTER(LEN=80)                :: chrkey, cprnt
  CHARACTER(LEN=3)                 :: itoc
  INTEGER                          :: intkey, chrlen, len
  INTEGER                          :: i, j, nply, nsubr, nblock, begin, ntest
  LOGICAL                          :: dollar, logkey, automte, mprnt
  REAL*8                           :: fpkey, step, boundl, boundr
  WRITE(output,1) coord
  LEN=chrlen(coord)
!
!     read the input by stopping at the keyword in the input stream.
!
  keywrd='$h0('//coord(1:LEN)//')'
  IF ( dollar(keywrd,card,cpass,input) ) THEN
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
!     or space coordinate
!
      IF(coord(1:1) /= 't') THEN
!
!       this is a space calculation
!
        WRITE(output,2)
        mass=fpkey(card,'mass',massau,' ')
        units=chrkey(card,'units','atomic-units',' ')
        parity=chrkey(card,'parity','none',' ')
        angmom=intkey(card,'angular-momentum',0,' ')
        reuse=logkey(card,'reuse-space-data',.false.,' ')
        IF(units == 'atomic-units') THEN
           mass=mass/massau
        END IF
        nodiag=logkey(card,'do-not-diagonalize',.false.,' ')
        automte=logkey(card,'automate',.false.,' ')
        typwt=chrkey(card,'weight-type','one',' ')
        IF(typwt == 'fourier') then
            nreg=1
            nfix=0
            fix(1)=.false.
            fix(2)=.false.
            drop(1)=.false.
            drop(2)=.false.
            bcl=1
            bcr=1
            edge(1)=-pi
            edge(2)=pi
            CALL fparr(card,'region-boundaries',edge,nreg+1,' ')
            CALL intarr(card,'number-of-points',n,nreg,' ')
            WRITE(output,3) edge(1), edge(2)
!
!           if n is not odd fix it.
!
            ntest = n(1) - 2 * (n(1)/2)
            IF(ntest == 0 ) THEN
               n(1) = n(1) + 1
            END IF 
            npt(1)=n(1)
            nrq(1)=npt(1) 
            WRITE(output,4) (npt(i),i=1,nreg)
        ELSE IF (typwt == 'legendre') then
            nreg=1
            m_val=intkey(card,'legendre-m',0,' ')
            nfix=0
            fix(1)=.false.
            fix(2)=.false.
            drop(1)=.false.
            drop(2)=.false.
            bcl=1
            bcr=1
            edge(1)=-1.d0
            edge(2)=1.d0
            IF(m_val /= 0) THEN
               nfix=2
               fix(1)=.true.
               fix(2)=.true.
               drop(1)=.true.
               drop(2)=.true.
               bcl=0
               bcr=0
            END IF 
            CALL intarr(card,'number-of-points',n,nreg,' ')
            npt(1)=n(1)
            nrq(1)=npt(1) 
            WRITE(output,4) (npt(i),i=1,nreg)
            write(output,5) m_val
        ELSE
            IF (typwt == 'hermite') then
                WRITE(output,6)
            END IF
            IF(typwt == 'laguerre') then
               WRITE(output,7)
            END IF
            nfix=intkey(card,'number-of-fixed-points',2,' ')
            write(output,8) nfix
            IF(nfix /= 0) THEN
               fix(1)=logkey(card,'left-fixed-point',.true.,' ')
               fix(2)=logkey(card,'right-fixed-point',.true.,' ')
               drop(1)=logkey(card,'drop-left-function',.false.,' ')
               drop(2)=logkey(card,'drop-right-function',.false.,' ')
               write(output,9) fix(1), drop(1), fix(2), drop(2)
            END IF
            bcl=1
            bcr=1
            IF(drop(1)) THEN
                bcl=0
            END IF
            IF(drop(2)) THEN
               bcr=0
            END IF
            IF(.not.automte) then
              write(output,10) 
              nreg=intkey(card,'number-of-regions',1,' ')
              CALL fparr(card,'region-boundaries',edge,nreg+1,' ')
              CALL intarr(card,'polynomial-order-per-region',n,nreg,' ')
              npt=n+1
              nrq=npt 
              CALL intarr(card,'number-of-reference-quadrature-'//  &
                               'points-per-region',nrq,nreg,' ')
              write(output,11) nreg
              write(output,4) (n(i),i=1,nreg)
              write(output,12) (edge(i),i=1,nreg+1)
              write(output,13) (nrq(i),i=1,nreg)
            ELSE
              nreg=0
              write(output,14)
              nblock=intkey(card,'number-of-major-blocks',1,' ')       
              DO i=1,nblock 
                 IF ( dollar('$block-'//itoc(i),card, &
                      cpass,input) ) THEN
                      skip=logkey(card,'skip',.false.,' ')
                      IF(skip) THEN
                         write(output,15) i
                      ELSE 
                         nply=intkey(card,'default-order',10,' ')
                         nsubr=intkey(card,'number-of-subregions',1,' ')
                         boundl=fpkey(card,'left-boundary',0.d0,' ')
                         boundr=fpkey(card,'right-boundary',0.d0,' ')
                         write(output,16) i, nply, nsubr, boundl, boundr
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
                 END IF
              END DO
              npt = n + 1
              nrq = npt 
            END IF
        END IF
!
!       Determine nphy
!
        call ptcal(nphy,nglobal,typwt)
        write(output,17) nphy
!
      ELSE IF (coord(1:1) == 't') THEN
!
!        this is a time calculation
!    
        WRITE(output,18)
        CALL fparr(card,'region-boundaries',edge,2,' ')
        endpts(1)=edge(1)
        endpts(2)=edge(2)
        n(1)=intkey(card,'polynomial-order',1,' ')
        typwt='one'
        reuse=logkey(card,'reuse-time-data',.false.,' ')
        WRITE(output,4) (n(i),i=1,2)
        WRITE(output,12) (edge(i),i=1,2)
!    
!       number of quadrature points is one greater
!       than polynomial order.
    
        npt(1)=n(1)+1
        nphy=n(1)
      ELSE
        write(output,19)
        stop 
      END IF
  END IF
1    FORMAT(/,20X,'coordinate = ',a24)
2    FORMAT(/,20X,'generating dvr representation in space')
3    FORMAT(/,1x,'     Fourier : Region is (',e15.8,',',e15.8,')')
4    FORMAT(/,15x,'polynomial order                       = ',  &
                   (/,15x,5(i4,1x)))                            
5    Format(/,1x,'legendre m value = ',i2)
6    FORMAT(/,1x,'     Hermite : Region is - infinity to + infinity')
7    FORMAT(/,1x,'     Laguerre: Region is zero to + infinity')
8    Format(/,1x,'number of fixed points = ',i1)
9    Format(/,15x,'left point fixed    = ',l1,/,15x,  &
                  'left point dropped  = ',l1,/,15x,  &
                  'right point fixed   = ',l1,/,15x,  &
                  'right point dropped = ',l1)
10   Format(/,1x,'point generation not automated')
11   FORMAT(/,1X,'number of regions = ', i3)
12   FORMAT(/,15x,'region boundaries                      = ',  &
                   (/,15x,5(e15.8,1x)))                         
13   FORMAT(/,15x,'number of reference quadrature points  = ',  &
                   (/,15x,5(i4,1x)))
14   Format(/,1x,'point generation automated')
15   FORMAT(/,1x,'skipping input block = ',i4)
16   FORMAT(/,1x, 'block = ',i3, &
            /,15x,'polynomial order     = ',i4, &
            /,15x,'number of subregions = ',i4, & 
            /,15x,'left hand boundary   = ',e15.8, &
            /,15x,'right hand boundary  = ',e15.8)
17   FORMAT(/,1x,'     size of the physical basis set = ',i5)
18   FORMAT(/,20X,'generating dvr representation in time')
19   FORMAT(/,1x,'no input found corresponding to keyword = ',a80)

END SUBROUTINE dvr_input

























