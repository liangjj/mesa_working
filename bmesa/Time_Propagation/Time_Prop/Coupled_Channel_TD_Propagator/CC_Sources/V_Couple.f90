!deck v_couple.f
!**begin prologue     v_couple
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            time-dependent potential
!**
!**description        input data for coupling potential, including the
!**                   mean field interaction for a contact potential.
!**references
!**routines called
!**end prologue       v_couple
  SUBROUTINE v_couple
  USE arnoldi_global
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                      :: natom, intkey
  REAL*8                       :: fpkey, sctlen
  CHARACTER (LEN=80)           :: chrkey
  CHARACTER (LEN=2)            :: atom
  LOGICAL                      :: dollar, fretyp, logkey
  IF(dollar('$v_couple',card,cpass,inp) ) THEN
     IF(system == 'spherical') THEN
        WRITE(iout,1)
        WRITE(iout,3)
        CALL intarr(card,'channel-l-values',l_chan,nchan,' ')
        CALL intarr(card,'channel-m-values',m_chan,nchan,' ')
        DO i=1,nchan
           WRITE(iout,4) l_chan(i), m_chan(i)
        END DO
        
     ELSE IF(system == 'cylindrical') THEN
        WRITE(iout,2)
        WRITE(iout,5)
        CALL intarr(card,'channel-m-values',m_chan,nchan,' ')
        DO i=1,nchan 
           WRITE(iout,6) m_chan(i)
        END DO
     END IF
  END IF
1    FORMAT(/,1X,'Spherical Coordinates')
2    FORMAT(/,1X,'Cylindrical Coordinates')
3    FORMAT(/,10x,'Channel L Value',10x,'Channel M Value')
4    FORMAT(16x,i3,23x,i3)
5    FORMAT(/,10x,'Channel M Value')
6    FORMAT(16x,i3)
END SUBROUTINE v_couple


