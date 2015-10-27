! \documentclass{article}
 
! Code converted using TO_F90 by Alan Miller
! Date: 2004-03-13  Time: 16:05:05
 
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{VTRSEP: Compute Davidson Trial Vectors}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle

!deck vtrsep.f
!***begin prologue     vtrsep
!***date written       010828   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           hamiltonian
!***author             schneider, barry (nsf)
!***source             3-dim
!***purpose            guess vectors based on separable hamiltonian
!***                   for davidson routine.
!***
!***
!***references

!***routines called
!***end prologue       vtrsep

SUBROUTINE vtrsep(u01,u02,u03,eig01,eig02,eig03,  &
    vec,eig,ind,n1,n2,n3,dim,n,nroots,prnt)

REAL*8, INTENT(IN)                       :: u01(n1,n1)
REAL*8, INTENT(IN)                       :: u02(n2,n2)
REAL*8, INTENT(IN OUT)                   :: u03(n3,n3)
REAL*8, INTENT(IN)                       :: eig01(n1)
REAL*8, INTENT(IN)                       :: eig02(n2)
REAL*8, INTENT(IN)                       :: eig03(n3)
REAL*8, INTENT(OUT)                      :: vec(n,nroots)
REAL*8, INTENT(OUT)                      :: eig(nroots)
INTEGER, INTENT(IN OUT)                  :: ind(n,*)
INTEGER, INTENT(IN)                      :: n1
INTEGER, INTENT(IN)                      :: n2
INTEGER, INTENT(IN)                      :: n3
INTEGER, INTENT(IN)                      :: dim
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN)                      :: nroots
LOGICAL, INTENT(IN)                      :: prnt
IMPLICIT INTEGER (a-z)


REAL*8  tmp, temp
CHARACTER (LEN=80) :: title, cpass, chrkey
CHARACTER (LEN=240) :: card
CHARACTER (LEN=5) :: itoc
CHARACTER (LEN=8) :: code
LOGICAL :: dollar



COMMON/io/inp, iout
DATA code / 'trial:' /

IF( dollar('$trial',card,cpass,inp) ) THEN
  title=chrkey(card,'type-trials','unperturbed',' ')
END IF

IF(title == 'unit') THEN
  CALL rzero(vec,n*nroots)
  DO  i=1,nroots
    vec(i,i)=1.d0
  END DO
ELSE IF(title == 'unperturbed') THEN
  
!       the ind array is destroyed by this routine and needs to be copied or
!       regenerated if it is used again.
  
  IF(dim == 1) THEN
    CALL copy(eig01,eig,nroots)
    CALL copy(u01,vec,n*nroots)
  ELSE IF(dim == 2) THEN
    count=0
    DO  i=1,n1
      DO  j=1,n2
        count=count+1
        eig(count) = eig01(i) + eig02(j)
      END DO
    END DO
    DO  ii=2,n
      i=ii-1
      k=i
      tmp=eig(i)
      i1=ind(i,1)
      j1=ind(i,2)
      DO  j=ii,n
        IF(eig(j) < tmp) THEN
          k=j
          tmp=eig(j)
        END IF
      END DO
      IF(k /= i) THEN
        ind(i,1)=ind(k,1)
        ind(i,2)=ind(k,2)
        ind(k,1)=i1
        ind(k,2)=j1
        eig(k) = eig(i)
        eig(i) = tmp
      END IF
    END DO
    DO  i=1,nroots
      count=0
      DO  j=1,n1
        DO  k=1,n2
          count=count+1
          vec(count,i)=u01(j,ind(i,1))*u02(k,ind(i,2))
        END DO
      END DO
    END DO
  ELSE
    count=0
    DO  i=1,n1
      DO  j=1,n2
        DO  k=1,n3
          count=count+1
          eig(count) = eig01(i) + eig02(j) + eig03(k)
        END DO
      END DO
    END DO
    DO  ii=2,n
      i=ii-1
      k=i
      tmp=eig(i)
      i1=ind(i,1)
      j1=ind(i,2)
      k1=ind(i,3)
      DO  j=ii,n
        IF(eig(j) < tmp) THEN
          k=j
          tmp=eig(j)
        END IF
      END DO
      IF(k /= i) THEN
        ind(i,1)=ind(k,1)
        ind(i,2)=ind(k,2)
        ind(i,3)=ind(k,3)
        ind(k,1)=i1
        ind(k,2)=j1
        ind(k,3)=k1
        eig(k) = eig(i)
        eig(i) = tmp
      END IF
    END DO
    DO  i=1,nroots
      count=0
      DO  j=1,n1
        DO  k=1,n2
          DO  l=1,n3
            count=count+1
            vec(count,i) = u01(j,ind(i,1)) * u02(k,ind(i,2)) *  &
                u03(l,ind(i,3))
          END DO
        END DO
      END DO
    END DO
  END IF
  IF(prnt) THEN
    title='guess eigenvalues'
    CALL prntfm(title,eig,nroots,1,nroots,1,iout)
  END IF
ELSE
  CALL lnkerr('error in trial type')
END IF
DO  i=1,nroots
  CALL iosys('write real "'//code(1:6)//itoc(i)//'" to rwf', n,vec(1,i),0,' ')
END DO
IF(prnt) THEN
  title='guess eigenvectors'
  CALL prntfm(title,vec,n,nroots,n,nroots,iout)
END IF
RETURN
END SUBROUTINE vtrsep

















