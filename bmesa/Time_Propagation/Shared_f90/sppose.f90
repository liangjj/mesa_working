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
  USE arnoldi_global
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  REAL*8, DIMENSION(:), ALLOCATABLE          :: cx, phix, cy, phiy, & 
                                                cz, phiz
  INTEGER, DIMENSION(:), ALLOCATABLE         :: lx, ly, lz
  LOGICAL                                    :: dollar
  INTEGER                                    :: ntot, nxs, nys, nzs, nrs
  INTEGER                                    :: intkey
  psi0=(0.d0,0.d0)
  IF ( dollar('$states',card,title,inp) ) THEN
       IF(system == 'cartesian') THEN
           ntot=1
           IF(spdim >= 1) THEN
              nxs=intkey(card,'number-of-x-superposed-states',1,' ')
              ntot=ntot*nxs
              ALLOCATE(cx(nxs),phix(nphy(1)),lx(nxs))
              CALL intarr(card,'x-state-list',lx,nxs,' ')
              CALL fparr(card,'x-state-coefficients',cx,nxs,' ')
              CALL mk_phi(phix,grid(1)%eigvec_0,cx,lx,nxs,nphy(1),'x')
           END IF
           IF(spdim >= 2) THEN
              nys=intkey(card,'number-of-y-superposed-states',1,' ')
              ntot=ntot*nys
              ALLOCATE(cy(nys),phiy(nphy(2)),ly(nys))
              CALL intarr(card,'y-state-list',ly,nys,' ')
              CALL fparr(card,'y-state-coefficients',cy,nys,' ')
              CALL mk_phi(phiy,grid(2)%eigvec_0,cy,ly,nys,nphy(2),'y')
           END IF
           IF(spdim == 3) THEN
              nzs=intkey(card,'number-of-z-superposed-states',1,' ')
              ntot=ntot*nzs
              ALLOCATE(cz(nzs),phiz(nphy(3)),lz(nzs))
              CALL intarr(card,'z-state-list',lz,nzs,' ')
              CALL fparr(card,'z-state-coefficients',cz,nzs,' ')
              CALL mk_phi(phiz,grid(3)%eigvec_0,cz,lz,nzs,nphy(3),'z')
           END IF
       ELSE
           nrs=intkey(card,'number-of-r-superposed-states',1,' ')
           ntot=nrs
           ALLOCATE(cx(nrs),phix(nphy(1)),lx(nrs))
           CALL intarr(card,'r-state-list',lx,nrs,' ')
           CALL fparr(card,'r-state-coefficients',cx,nrs,' ')
           CALL mk_phi(phix,grid(1)%eigvec_0,cx,lx,nrs,nphy(1),'r')
       END IF
  END IF
  IF(spdim == 1) THEN
     call fil_1(phix,nphy(1))
     DEALLOCATE(cx,phix,lx)
  END IF
  IF(spdim == 2) THEN
     CALL fil_2(phix,phiy,nphy(1),nphy(2))
     DEALLOCATE(cx,phix,lx)
     DEALLOCATE(cy,phiy,ly)
  END IF
  IF(spdim == 3) THEN
     CALL fil_3(phix,phiy,phiz,nphy(1),nphy(2),nphy(3))
     DEALLOCATE(cx,phix,lx)
     DEALLOCATE(cy,phiy,ly)
     DEALLOCATE(cz,phiz,lz)
  END IF
END SUBROUTINE sppose
