! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Generalized Kinetic Energy Matrix Elements}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck kemel.f
!***begin prologue     kemel
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           coordinate eigenfunctions
!***author             schneider, barry (nsf)
!***source
!***purpose            generate kinetic energy matrix elements with singularities
!***                   removed and bloch operators added.
!***
!***description
!***references
!***routines called
!***end prologue       kemel
  SUBROUTINE kemel(kmat,p,dp,ddp,pt,wt,edge,n, cordsys,coord,parity)
  USE dvr_global,    ONLY   : inp,iout
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n,n)                 :: kmat, p, dp, ddp
  REAL*8, DIMENSION(n)                   :: pt, wt
  REAL*8, DIMENSION(2)                   :: edge
  CHARACTER (LEN=*)                      :: cordsys, coord, parity
  CHARACTER (LEN=80)                     :: title
  kmat=0.d0
  IF(cordsys == 'cartesian') THEN
     IF(coord == 'x'.OR. coord == 'y'.OR.  &
        coord == 'z') THEN
! The kinetic energy operator is assumed to be,
! \begin{equation}
!         T = - \frac{1}{2M} \frac{d^2}{dx^2}
! \end{equation}
! $x$ can be the transformed radial equation as well as cartesian coordinates.
       CALL kinxyz(kmat,p,ddp,wt,n)
       CALL blxyz(kmat,p,dp,pt(1),edge(1),pt(n),edge(2),n)
     ELSE IF(coord == 'r') THEN
       CALL kinxyz(kmat,p,ddp,wt,n)
       CALL blxyz(kmat,p,dp,pt(1),edge(1),pt(n),edge(2),n)
     ELSE
       CALL lnkerr('error in axis type')
     END IF
  ELSE IF(cordsys == 'spherical') THEN
! The kinetic energy operator is assumed to be,
! \begin{equation}
!    T =- \frac{1}{2M} \frac{1}{r^2} \frac{d}{dr} \big (  r^2 \frac{d} {dr} \big )
! \end{equation}
     CALL kinrad(kmat,p,dp,ddp,pt,wt,n,parity)
     CALL blrad(kmat,p,dp,pt(1),edge(1),pt(n),edge(2),n,parity)
  ELSE IF(cordsys == 'cylindrical') THEN
! The kinetic energy operator is assumed to be,
! \begin{equation}
!    T =- \frac{1}{2M} \frac{1}{\rho} \frac{d}{d \rho} \big (  \rho \frac{d} {d \rho} \big )
! \end{equation}
     CALL kincyl(kmat,p,dp,ddp,pt,wt,n,parity)
     CALL blcyl(kmat,p,dp,pt(1),edge(1),pt(n),edge(2),n,parity)
  ELSE
     CALL lnkerr('quit')
  END IF
END SUBROUTINE kemel
