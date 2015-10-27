! \documentclass{article}
 
! Code converted using TO_F90 by Alan Miller
! Date: 2008-03-16  Time: 12:02:45
 
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{GENVEC: Generate New Trial Vectors based on Diagonal Preconditioning}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle

!deck genvec.f
!***begin prologue     genvec
!***date written       010829   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           davidson, trial, vector
!***author             schneider, barry (nsf)
!***source
!***purpose            generate a new trial vector based on some zeroth
!***                   order model.  this is equivalent to a preconditioning
!***                   of the matrix.
!***references

!***routines called
!***end prologue       genvec

SUBROUTINE genvec(vec,resid,diag,eig,work,n1,n2,n3,n,dim,  &
    p,add,iter,precon,maxblk,prnt)

REAL*8, INTENT(OUT)                      :: vec(n,add)
REAL*8, INTENT(IN)                       :: resid(n,add)
REAL*8, INTENT(IN)                       :: diag(n)
REAL*8, INTENT(IN)                       :: eig(add)
REAL*8, INTENT(IN OUT)                   :: work(*)
INTEGER, INTENT(IN)                      :: n1
INTEGER, INTENT(IN)                      :: n2
INTEGER, INTENT(IN)                      :: n3
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN)                      :: dim
INTEGER*8, INTENT(IN OUT)                :: p
INTEGER, INTENT(IN)                      :: add
INTEGER, INTENT(IN OUT)                  :: iter
CHARACTER (LEN=*), INTENT(IN)            :: precon
INTEGER, INTENT(IN)                      :: maxblk
LOGICAL, INTENT(IN)                      :: prnt
IMPLICIT INTEGER (a-z)

REAL*8 test, zero, nrzero, one, array
#ifdef decpointer

#END IF decpointer
#ifdef sgipointer

#END IF sgipointer

CHARACTER (LEN=4) :: itoc
CHARACTER (LEN=80) :: title

DATA zero, nrzero, one / 0.d0, 1.0D-06, 1.d0 /

COMMON/io/inp, iout
pointer (p,array(1)), (p,iarray(1))

IF(prnt) THEN
  title='input residuals to genvec'
  CALL prntrm(title,resid,n,add,n,add,iout)
END IF
IF(precon == 'diagonal') THEN
  DO  i=1,add
    DO  j=1,n
      test=eig(i) - diag(j)
      IF(ABS(test) >= nrzero) THEN
        vec(j,i) = resid(j,i)/test
      ELSE
        vec(j,i) = one
      END IF
    END DO
  END DO
!         title='vectors'
!         call prntrm(title,vec,n1,add,n1,add,iout)
ELSE IF(precon == 'separable') THEN
  u1=1
  eig1=u1+n1*n1
  u2=eig1+n1
  eig2=u2+n2*n2
  u3=eig2+n2
  eig3=u3+n3*n3
  title='eigenvalues'
  CALL prntrm(title,array(eig1),n1,1,n1,1,iout)
  title='eigenvalues'
  CALL prntrm(title,eig,add,1,add,1,iout)
  title='transformation matrix'
  CALL prntrm(title,array(u1),n1,n1,n1,n1,iout)
  title='residuals'
  CALL prntrm(title,resid,n1,add,n1,add,iout)
  IF(dim == 1) THEN
    CALL h12h0(resid,work,array(u1),n1,add)
    title='h12h0 work'
    CALL prntrm(title,work,n1,add,n1,add,iout)
    CALL h1e(work,array(eig1),eig,n1,add)
    title='h1e work'
    CALL prntrm(title,work,n1,add,n1,add,iout)
    CALL h12dvr(work,vec,array(u1),n1,add)
    title='vectors'
    CALL prntrm(title,vec,n1,add,n1,add,iout)
  ELSE IF(dim == 2) THEN
    CALL h22h0(resid,work,array(u1),array(u2),n1,n2,add)
    CALL h2e(work,array(eig1),array(eig2),eig,n1,n2,add)
    CALL h22dvr(work,vec,array(u1),array(u2),n1,n2,add)
  ELSE IF(dim == 3) THEN
    CALL h32h0(resid,work,array(u1),array(u2),array(u3), n1,n2,n3,add)
    CALL h3e(work,array(eig1),array(eig2),array(eig3),eig, n1,n2,n3,,add)
    CALL h32dvr(work,vec,array(u1),array(u2),array(u3), n1,n2,n3,add)
  END IF
ELSE IF(precon == 'block') THEN
  u0=1
  eig0=u0+maxblk*maxblk
  CALL rblvec(resid,vec,array(eig0),eig,array(u0),n,add)
END IF
IF(prnt) THEN
  title='new trial vectors iteration = '//itoc(iter)
  CALL prntrm(title,vec,n,add,n,add,iout)
END IF
RETURN
END SUBROUTINE genvec





