!deck modify_diag.f
!***begin prologue     modify_diag
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            modify diagonal DVR matrix elements
!***
!***references
!***routines called
!***end prologue       modify_diag
!
  SUBROUTINE modify_diag(ke,v,nr,nf,modify,mat_typ)
  USE dvr_global
  USE dvrprop_global
  USE fd_global
  IMPLICIT NONE
  INTEGER                                :: nr, nf, i
  REAL*8, DIMENSION(nf)                  :: v
  REAL*8, DIMENSION(nr,nf)               :: ke
  CHARACTER(LEN=*)                       :: mat_typ, modify
!
!
  IF(modify == 'zero_potential_energy') THEN
!
!    here we can zero the potential energy if it is included
!    somewhere else in the problem.
!
     v=0.d0
     write(iout,1)
  ELSE
     IF( mat_typ == 'full') THEN
         IF(modify == 'add_potential_energy') THEN
!
!           add the potential energy to the kinetic energy.
!           then zero it.
!
            write(iout,2)
            DO i=1,nf
               ke(i,i) = ke(i,i) + v(i)
               v(i) = ke(i,i)
            END DO
            if(log_main(1)) then
                title='Modified Diagonal of Kinetic Energy Matrix. Potential Set'  &
                       //' To Zero'
                call prntfmn(title,v,nf,1,nf,1,iout,'e')
            END IF
            v = 0.d0
         ELSE IF(modify == 'zero_kinetic_energy') THEN
            write(iout,3)
!
!           absorb the diagonal kinetic energy into the potential
!           then zero the diagonal kinetic energy.
!
            DO i=1,nf
               v(i) = ke(i,i) + v(i)
            END DO
            if(log_main(1)) then
               title='Modified Potential Energy Matrix. Diagonal Kinetic Energy ' &
                     //'Set To Zero'
               call prntfmn(title,v,nf,1,nf,1,iout,'e')
           END IF
           DO i=1,nf
              ke(i,i) = 0.d0
           END DO
         END IF
     ELSE IF( mat_typ == 'banded') THEN
         IF(modify == 'add_potential_energy') THEN
            DO i=1,nf
               ke(1,i) = ke(1,i) + v(i)
               v(i) = ke(1,i)
            END DO
            if(log_main(1)) then
                title='Modified Diagonal of Kinetic Energy Matrix. Potential Set '  &
                       //'To Zero'
               call prntfmn(title,v,nf,1,nf,1,iout,'e')
            END IF
            v = 0.d0
         ELSE IF(modify == 'zero_kinetic_energy') THEN
            DO i=1,nf
               v(i) = v(i) + ke(1,i)
            END DO
            if(log_main(1)) then
                title='Modified Potential Energy Matrix. Diagonal Kinetic Energy '  &
                       //'Set To Zero'
                call prntfmn(title,v,nf,1,nf,1,iout,'e')
            END IF
            ke(1,:) = 0.d0
         END IF
     ELSE
        call lnkerr('quit')
     END IF
  END IF
1 FORMAT(/,10x,'V Zeroed and KE Left Unmodified')
2 FORMAT(/,10x,'Diagonal of KE Modified to Include One Body V and Then' &
               ' V Zeroed')
3 FORMAT(/,10x,'V Modified to Include Diagonal KE and Then Diagonal KE' &
               ' Zeroed')
END SUBROUTINE modify_diag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
