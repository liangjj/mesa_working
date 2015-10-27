! \documentclass{article}
 
! Code converted using TO_F90 by Alan Miller
! Date: 2004-03-13  Time: 16:04:55
 
! \ublkackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{VTRBLK: Compute Davidson Trial Vectors}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle

!deck vtrblk.f
!***begin prologue     vtrblk
!***date written       010828   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           hamiltonian
!***author             schneider, barry (nsf)
!***source             3-dim
!***purpose            guess vectors based on blocked hamiltonian
!***                   for davidson routine.
!***
!***
!***references

!***routines called
!***end prologue       vtrblk

SUBROUTINE vtrblk(eig,eigtr,vec,ind,n,ntrial,prnt)

REAL*8, INTENT(IN OUT)                   :: eig(n)
REAL*8, INTENT(IN OUT)                   :: eigtr(ntrial)
REAL*8, INTENT(IN OUT)                   :: vec(n,ntrial)
INTEGER, INTENT(OUT)                     :: ind(n)
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN)                      :: ntrial
LOGICAL, INTENT(IN)                      :: prnt
IMPLICIT INTEGER (a-z)


CHARACTER (LEN=24) :: str
CHARACTER (LEN=2) :: itoc
CHARACTER (LEN=8) :: code
CHARACTER (LEN=80) :: title


COMMON/io/inp, iout
DATA code / 'trial:' /

CALL iosys('rewind vectors on ham read-and-write',0,0,0,' ')
DO  i=1,n
  ind(i)=i
END DO
DO  ii=2,n
  i=ii-1
  k=i
  tmp=eig(i)
  i1=ind(i)
  DO  j=ii,n
    IF(eig(j) < tmp) THEN
      k=j
      tmp=eig(j)
    END IF
  END DO
  IF(k /= i) THEN
    ind(i)=ind(k)
    ind(k)=i1
    eig(k) = eig(i)
    eig(i) = tmp
  END IF
  CALL copy(eig,eigtr,ntrial)
  DO  trl=1,ntrial
    itr=ind(trl)
    offset=(itr-1)*n
    CALL iosys('read real vector from ham',n,vec(1,trl), offset,' ')
  END DO
END DO
IF(prnt) THEN
  title='guess eigenvalues'
  CALL prntfm(title,eig,nroots,1,nroots,1,iout)
END IF
DO  i=1,nroots
  CALL iosys('write real "'//code(1:6)//itoc(i)//'" to rwf', n,vec(1,i),0,' ')
END DO
IF(prnt) THEN
  title='guess eigenvectors'
  CALL prntfm(title,vec,n,nroots,n,nroots,iout)
END IF
RETURN
END SUBROUTINE vtrblk

















