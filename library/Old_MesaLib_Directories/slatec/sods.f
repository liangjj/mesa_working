*deck sods
      subroutine sods (a, x, b, neq, nuk, nrda, iflag, work, iwork)
c***begin prologue  sods
c***subsidiary
c***purpose  subsidiary to bvsup
c***library   slatec
c***type      single precision (sods-s)
c***author  watts, h. a., (snla)
c***description
c
c     sods solves the overdetermined system of linear equations a x = b,
c     where a is neq by nuk and neq .ge. nuk. if rank a = nuk,
c     x is the unique least squares solution vector. that is,
c              r(1)**2 + ..... + r(neq)**2 = minimum
c     where r is the residual vector  r = b - a x.
c     if rank a .lt. nuk , the least squares solution of minimal
c     length can be provided.
c     sods is an interfacing routine which calls subroutine lssods
c     for the solution. lssods in turn calls subroutine orthol and
c     possibly subroutine ohtror for the decomposition of a by
c     orthogonal transformations. in the process, orthol calls upon
c     subroutine cscale for scaling.
c
c **********************************************************************
c   input
c **********************************************************************
c
c     a -- contains the matrix of neq equations in nuk unknowns and must
c          be dimensioned nrda by nuk. the original a is destroyed
c     x -- solution array of length at least nuk
c     b -- given constant vector of length neq, b is destroyed
c     neq -- number of equations, neq greater or equal to 1
c     nuk -- number of columns in the matrix (which is also the number
c            of unknowns), nuk not larger than neq
c     nrda -- row dimension of a, nrda greater or equal to neq
c     iflag -- status indicator
c            =0 for the first call (and for each new problem defined by
c               a new matrix a) when the matrix data is treated as exact
c           =-k for the first call (and for each new problem defined by
c               a new matrix a) when the matrix data is assumed to be
c               accurate to about k digits
c            =1 for subsequent calls whenever the matrix a has already
c               been decomposed (problems with new vectors b but
c               same matrix a can be handled efficiently)
c     work(*),iwork(*) -- arrays for storage of internal information,
c                     work must be dimensioned at least  2 + 5*nuk
c                     iwork must be dimensioned at least nuk+2
c     iwork(2) -- scaling indicator
c                 =-1 if the matrix a is to be pre-scaled by
c                 columns when appropriate
c                 if the scaling indicator is not equal to -1
c                 no scaling will be attempted
c              for most problems scaling will probably not be necessary
c
c **********************************************************************
c   output
c **********************************************************************
c
c     iflag -- status indicator
c            =1 if solution was obtained
c            =2 if improper input is detected
c            =3 if rank of matrix is less than nuk
c               if the minimal length least squares solution is
c               desired, simply reset iflag=1 and call the code again
c     x -- least squares solution of  a x = b
c     a -- contains the strictly upper triangular part of the reduced
c           matrix and the transformation information
c     work(*),iwork(*) -- contains information needed on subsequent
c                         calls (iflag=1 case on input) which must not
c                         be altered
c                         work(1) contains the euclidean norm of
c                         the residual vector
c                         work(2) contains the euclidean norm of
c                         the solution vector
c                         iwork(1) contains the numerically determined
c                         rank of the matrix a
c
c **********************************************************************
c
c***see also  bvsup
c***references  g. golub, numerical methods for solving linear least
c                 squares problems, numerische mathematik 7, (1965),
c                 pp. 206-216.
c               p. businger and g. golub, linear least squares
c                 solutions by householder transformations, numerische
c                 mathematik  7, (1965), pp. 269-276.
c               h. a. watts, solving linear least squares problems
c                 using sods/suds/cods, sandia report sand77-0683,
c                 sandia laboratories, 1977.
c***routines called  lssods
c***revision history  (yymmdd)
c   750601  date written
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   910408  updated the author and references sections.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  sods
      dimension a(nrda,*),x(*),b(*),work(*),iwork(*)
c
c***first executable statement  sods
      iter=0
      is=2
      ip=3
      ks=2
      kd=3
      kz=kd+nuk
      kv=kz+nuk
      kt=kv+nuk
      kc=kt+nuk
c
      call lssods(a,x,b,neq,nuk,nrda,iflag,iwork(1),iwork(is),a,
     1            work(kd),iwork(ip),iter,work(1),work(ks),
     2            work(kz),b,work(kv),work(kt),work(kc))
c
      return
      end
