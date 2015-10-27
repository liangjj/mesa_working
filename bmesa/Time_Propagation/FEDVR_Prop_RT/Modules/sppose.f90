!deck sppose.f
!***begin prologue     sppose
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose            calculate zero time gaussian wavepacket
!***
!***references
!***routines called
!***end prologue       sppose
  SUBROUTINE sppose
  USE dvrprop_global_rt
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  LOGICAL                                    :: dollar
  INTEGER                                    :: ntot
  INTEGER, DIMENSION(3)                      :: n_val
  INTEGER                                    :: i, intkey
  ALLOCATE(c_loc(spdim))
  DO i=1,spdim
     ALLOCATE(c_loc(i)%c(nphy(i)),c_loc(i)%phi(nphy(i)),c_loc(i)%l(nphy(i)))
  END DO
  psi = 0.d0
  IF ( dollar('$states',card,title,inp) ) THEN
       IF(system == 'cartesian') THEN
           ntot=1
           IF(spdim >= 1) THEN
              n_val(1)=intkey(card,'number-of-x-superposed-states',1,' ')
              ntot=ntot*n_val(1)
              CALL intarr(card,'x-state-list',c_loc(1)%l,n_val(1),' ')
              CALL fparr(card,'x-state-coefficients',c_loc(1)%c,n_val(1),' ')
              CALL mk_phi(c_loc(1)%phi,grid(1)%eigvec_0,c_loc(1)%c,c_loc(1)%l, &
                          n_val(1),nphy(1),'x')
           END IF
           IF(spdim >= 2) THEN
              n_val(2)=intkey(card,'number-of-y-superposed-states',1,' ')
              ntot=ntot*n_val(2)
              CALL intarr(card,'y-state-list',c_loc(2)%l,n_val(2),' ')
              CALL fparr(card,'y-state-coefficients',c_loc(2)%c,n_val(2),' ')
              CALL mk_phi(c_loc(2)%phi,grid(2)%eigvec_0,c_loc(2)%c,c_loc(2)%l, &
                          n_val(2),nphy(2),'y')
           END IF
           IF(spdim == 3) THEN
              n_val(3)=intkey(card,'number-of-z-superposed-states',1,' ')
              ntot=ntot*n_val(3)
              CALL intarr(card,'z-state-list',c_loc(3)%l,n_val(3),' ')
              CALL fparr(card,'z-state-coefficients',c_loc(3)%c,n_val(3),' ')
              CALL mk_phi(c_loc(3)%phi,grid(3)%eigvec_0,c_loc(3)%c,c_loc(3)%l, &
                          n_val(3),nphy(3),'z')
           END IF
       ELSE
           n_val(1)=intkey(card,'number-of-r-superposed-states',1,' ')
           ntot=n_val(1)
           CALL intarr(card,'r-state-list',c_loc(1)%l,n_val(1),' ')
           CALL fparr(card,'r-state-coefficients',c_loc(1)%c,n_val(1),' ')
           CALL mk_phi(c_loc(1)%phi,grid(1)%eigvec_0,c_loc(1)%c,c_loc(1)%l, &
                       n_val(1),nphy(1),'r')
       END IF
  END IF
  call fill_d(c_loc(1)%phi,c_loc(2)%phi,c_loc(3)%phi,nphy)
  DO i=1,spdim
     DEALLOCATE(c_loc(i)%c,c_loc(i)%phi,c_loc(i)%l)
  END DO
  DEALLOCATE(c_loc)
END SUBROUTINE sppose
