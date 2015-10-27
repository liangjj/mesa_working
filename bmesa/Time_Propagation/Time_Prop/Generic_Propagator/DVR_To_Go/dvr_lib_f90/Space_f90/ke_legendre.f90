!deck ke_legendre.f
!***begin prologue     ke_legendre
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           coordinate eigenfunctions
!***author             schneider, barry (nsf)
!***source
!***purpose            generate kinetic energy matrix elements for
!***                   legendre on (-1,1).
!***
!***description
!***references
!***routines called
!***end prologue       ke_legendre
!\begin{eqnarray}
!\end{eqnarray}
  SUBROUTINE ke_legendre(kin,f,df,ddf,pt,wt,n,region)
  USE dvr_global,     ONLY   : parity, mass, iout, m_val
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n,n)                 :: kin, f, df, ddf
  REAL*8, DIMENSION(n)                   :: pt, wt
  INTEGER                                :: region
  REAL*8                                 :: scale
  CHARACTER (LEN=80)                     :: title
  CHARACTER (LEN=3)                      :: itoc
  INTEGER                                :: i, j, k
!  Code using second derivative and Bloch operator
!
  kin=0.d0
  DO i=1,n
     kin(i,:) = kin(i,:) + wt(i) * f(i,i) *                   &
                         ( ( 1.d0-pt(i)*pt(i) ) * ddf(i,:)    &
                            - 2.d0 * pt(i) *df(i,:) )
  END DO
  kin(n,:) = kin(n,:) - f(n,n) * ( 1.d0 - pt(n)*pt(n) ) * df(n,:)
  kin(1,:) = kin(1,:) + f(1,1) * ( 1.d0 - pt(1)*pt(1) ) * df(1,:)
  IF(m_val /= 0) then
     DO i=1,n
        kin(i,i) = kin(i,i) - wt(i) * f(i,i) * f(i,i) *     &
                              (m_val*m_val/(1.d0 - pt(i)*pt(i)))
     END DO
  END IF
  scale= -.5D0/mass
  kin=scale*kin
  IF(prn(3)) THEN
     title='unnormalized kinetic energy matrix for '// 'region = '//itoc(region)
     CALL prntrm(title,kin,n,n,n,n,iout)
  END IF
!  Code using first derivatives and no Bloch operator
!
!  kin=0.d0
!  DO i=1,n
!     DO j=1,n
!        DO k=1,n
!           kin(i,j) = kin(i,j) - df(k,i) * wt(k) *            &
!                         ( 1.d0-pt(k)*pt(k) ) * df(k,j)
!        END DO
!     END DO
!  END DO
!  IF(m_val /= 0) then
!     DO i=1,n
!        kin(i,i) = kin(i,i) - wt(i) * f(i,i) * f(i,i) *     &
!                              (m_val*m_val/(1.d0 - pt(i)*pt(i)))
!     END DO
!  END IF
!  scale= -.5D0/mass
!  kin=scale*kin
!  IF(prn(3)) THEN
!     title='unnormalized kinetic energy matrix for '// 'region = '//itoc(region)
!     CALL prntrm(title,kin,n,n,n,n,iout)
!  END IF
END SUBROUTINE ke_legendre
