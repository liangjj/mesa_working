!deck lares.f
!***begin prologue     lares
!***date written       980420   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           residual calculation
!***author             schneider, barry (nsf)
!***source
!***
!***references

!***routines called
!***end prologue       lares

SUBROUTINE lares(vec,hvec,b,coef,rhs,energy,cnverg,resid,  &
    maxerr,t,list,n,m,nrhs,ncon,addvec,maxvec, it,prnt)

REAL*8, INTENT(IN OUT)                   :: vec(n,maxvec)
REAL*8, INTENT(IN OUT)                   :: hvec(n,maxvec)
REAL*8, INTENT(IN OUT)                   :: b(maxvec,nrhs)
REAL*8, INTENT(IN OUT)                   :: coef(maxvec,nrhs)
REAL*8, INTENT(IN)                       :: rhs(n,nrhs)
REAL*8, INTENT(IN)                       :: energy
REAL*8, INTENT(IN)                       :: cnverg
REAL*8, INTENT(OUT)                      :: resid(n,nrhs)
REAL*8, INTENT(OUT)                      :: maxerr
REAL*8, INTENT(IN)                       :: t(n,nrhs)
INTEGER, INTENT(OUT)                     :: list(*)
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN OUT)                  :: m
INTEGER, INTENT(IN)                      :: nrhs
INTEGER, INTENT(OUT)                     :: ncon
INTEGER, INTENT(OUT)                     :: addvec
INTEGER, INTENT(IN OUT)                  :: maxvec
INTEGER, INTENT(IN OUT)                  :: it
LOGICAL, INTENT(IN)                      :: prnt(4)
IMPLICIT INTEGER (a-z)

REAL*8  sdot, ERR
CHARACTER (LEN=16) :: STATUS
CHARACTER (LEN=80) :: title
CHARACTER (LEN=4) :: itoc




COMMON/io/inp, iout

!        first express the solution vectors, Psi, in the original basis
!        and put them in work_d.
!
  CALL ebcxx(work_d,vec_d,small_rhs_work_d,matrix_size,end,1,                &
             matrix_size,matrix_size,maximum_number_of_davidson_vectors)
!
!        calculate H*Psi and put in residual_d.
!
  CALL ebcxx(residual_d,h_vectors_d,small_rhs_work_d,matrix_size,end,1,      &
             matrix_size,matrix_size,maximum_number_of_davidson_vectors)
!
  IF(prnt(1)) THEN
     title='information for iteration = '//itoc(it)
     WRITE(iout,1) title
  END IF
!
!        final residual
!
  residual_d(1:matrix_size) = rhs_d(1_matrix_size) - residual_d(1:matrix_size)
  IF(prnt(2)) THEN
     title='residuals iteration = '//itoc(it)
     CALL prntfm(title,residual_d,matrix_size,1,matrix_size,                 &
                 maximum_number_of_davidson_vectors,iout)
  END IF
  addvec=0
  ncon=0
  maxerr=0.d0
  ERR = SQRT (ddot(matrix_size,residual_d,1,residual_d,1) )
  maxerr=MAX(ERR,maxerr)

  IF(ERR <= cnverg) THEN
    STATUS='converged'
  ELSE
    STATUS='unconverged'    
    addvec = 1
  END IF
  WRITE(iout,2) ERR, STATUS
1    FORMAT(/,5X,a80)
2    FORMAT(/,10x,'RMS Error       = ',e15.8,/,5X, 'Status          = ',a16)
END SUBROUTINE lares
