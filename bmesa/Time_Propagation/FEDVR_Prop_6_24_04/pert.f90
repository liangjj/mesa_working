! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Pert}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck pert.f
!**begin prologue     pert
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            time-dependent potential
!**description        calculate the time and space dependent potential
!**                   matrix elements in the dvr representation.
!**references
!**routines called
!**end prologue       pert
  SUBROUTINE pert
  USE arnoldi_global
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  REAL*8                                 :: tfac
  INTEGER                                :: i, j, k, count
  CHARACTER(LEN=16)                      :: fptoc
  WRITE(iout,1)
!
!     There could be one body potentials.  Put them in v_tot.
  IF(spdim == 1 ) then
     v_tot = grid(1)%v
!     title='points'
!     call prntrm(title,grid(1)%pt,n3d,1,n3d,1,iout)
!     title='atomic potential'
!     call prntrm(title,v_tot,n3d,1,n3d,1,iout)
  ELSE IF(spdim == 2 ) then
     count = 0
     do i=1,nphy(1)
        do j=1,nphy(2)
           count = count + 1
           v_tot(count) = grid(1)%v(i) + grid(2)%v(j)
        END DO
     END DO
  ELSE IF(spdim == 3 ) then
     count = 0
     do i=1,nphy(1)
        do j=1,nphy(2)
           do k=1,nphy(3)
              count = count + 1
              v_tot(count) = grid(1)%v(i) + grid(2)%v(j) + grid(3)%v(k)
           END DO
        END DO
     END DO
  END IF
!     Deal with pure time potentials first.  They just need to be added it as
!     they were computed in a previous call.

  v_tot = v_tot + vt
!
!     The current routine has a number of options for potentials.  They include
!     sine, cosine and dipole-like terms.  Others could be added by the user.
!
  IF(vtyp(1) == 'cosine') THEN
      tfac = scale * COS(omega*t0)
      v_tot = v_tot + tfac
  ELSE IF(vtyp(1) == 'sine') THEN
      tfac = scale * SIN(omega*t0)
      v_tot = v_tot + tfac
  ELSE IF(vtyp(1) == 't') THEN
      v_tot = v_tot + t0
  ELSE IF(vtyp(1) == 'cosine-dipole-field') THEN
      tfac = scale * COS(omega*t0)
      CALL sc_dip(tfac)
  ELSE IF(vtyp(1) == 'sine-dipole-field') THEN
      tfac = scale * SIN(omega*t0)
      CALL sc_dip(tfac)
  ELSE IF(vtyp(1) == 'cosine-gaussian-pulse') THEN
      tfac = scale * COS(omega*t0)*EXP(-width*(t0-shift)*(t0-shift))
      CALL sc_dip(tfac)
  ELSE IF(vtyp(1) == 'sine-gaussian-pulse') THEN
      tfac = scale * COS(omega*t0) * EXP(-width*(t0-shift)*(t0-shift))
      CALL sc_dip(tfac)
  END IF
  IF(log_main(3)) THEN
     title='time-dependent perturbation for t = '//fptoc(t0) &
            //' omega = '//fptoc(omega)
     CALL prntfmn(title,v_tot,n3d,1,n3d,1,iout,'e')
  END IF
1    FORMAT(/,'calculating linear perturbing potential')
END SUBROUTINE pert


























