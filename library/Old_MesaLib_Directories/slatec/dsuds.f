*deck dsuds
      subroutine dsuds (a, x, b, neq, nuk, nrda, iflag, mlso, work,
     +   iwork)
c***begin prologue  dsuds
c***subsidiary
c***purpose  subsidiary to dbvsup
c***library   slatec
c***type      double precision (suds-s, dsuds-d)
c***author  watts, h. a., (snla)
c***description
c
c     dsuds solves the underdetermined system of linear equations a z =
c     b where a is neq by nuk and neq .le. nuk. in particular, if rank
c     a equals ira, a vector x and a matrix u are determined such that
c     x is the unique solution of smallest length, satisfying a x = b,
c     and the columns of u form an orthonormal basis for the null
c     space of a, satisfying a u = 0 .  then all solutions z are
c     given by
c              z = x + c(1)*u(1) + ..... + c(nuk-ira)*u(nuk-ira)
c     where u(j) represents the j-th column of u and the c(j) are
c     arbitrary constants.
c     if the system of equations are not compatible, only the least
c     squares solution of minimal length is computed.
c     dsuds is an interfacing routine which calls subroutine dlssud
c     for the solution.  dlssud in turn calls subroutine dorthr and
c     possibly subroutine dohtrl for the decomposition of a by
c     orthogonal transformations.  in the process, dorthr calls upon
c     subroutine dcscal for scaling.
c
c ********************************************************************
c   input
c ********************************************************************
c
c     a -- contains the matrix of neq equations in nuk unknowns and must
c          be dimensioned nrda by nuk.  the original a is destroyed.
c     x -- solution array of length at least nuk.
c     b -- given constant vector of length neq, b is destroyed.
c     neq -- number of equations, neq greater or equal to 1.
c     nuk -- number of columns in the matrix (which is also the number
c            of unknowns), nuk not smaller than neq.
c     nrda -- row dimension of a, nrda greater or equal to neq.
c     iflag -- status indicator
c           =0  for the first call (and for each new problem defined by
c               a new matrix a) when the matrix data is treated as exact
c           =-k for the first call (and for each new problem defined by
c               a new matrix a) when the matrix data is assumed to be
c               accurate to about k digits.
c           =1  for subsequent calls whenever the matrix a has already
c               been decomposed (problems with new vectors b but
c               same matrix a can be handled efficiently).
c     mlso -- =0 if only the minimal length solution is wanted.
c             =1 if the complete solution is wanted, includes the
c                linear space defined by the matrix u in the abstract.
c     work(*),iwork(*) -- arrays for storage of internal information,
c                work must be dimensioned at least
c                       nuk + 3*neq + mlso*nuk*(nuk-rank a)
c                where it is possible for   0 .le. rank a .le. neq
c                iwork must be dimensioned at least   3 + neq
c     iwork(2) -- scaling indicator
c                 =-1 if the matrix is to be pre-scaled by
c                 columns when appropriate.
c                 if the scaling indicator is not equal to -1
c                 no scaling will be attempted.
c              for most problems scaling will probably not be necessary
c
c *********************************************************************
c   output
c *********************************************************************
c
c     iflag -- status indicator
c            =1 if solution was obtained.
c            =2 if improper input is detected.
c            =3 if rank of matrix is less than neq.
c               to continue simply reset iflag=1 and call dsuds again.
c            =4 if the system of equations appears to be inconsistent.
c               however, the least squares solution of minimal length
c               was obtained.
c     x -- minimal length least squares solution of  a x = b.
c     a -- contains the strictly upper triangular part of the reduced
c           matrix and transformation information.
c     work(*),iwork(*) -- contains information needed on subsequent
c                         calls (iflag=1 case on input) which must not
c                         be altered.
c                         the matrix u described in the abstract is
c                         stored in the  nuk*(nuk-rank a) elements of
c                         the work array beginning at work(1+nuk+3*neq).
c                         however u is not defined when mlso=0 or
c                         iflag=4.
c                         iwork(1) contains the numerically determined
c                         rank of the matrix a
c
c *********************************************************************
c
c***see also  dbvsup
c***references  h. a. watts, solving linear least squares problems
c                 using sods/suds/cods, sandia report sand77-0683,
c                 sandia laboratories, 1977.
c***routines called  dlssud
c***revision history  (yymmdd)
c   750601  date written
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910408  updated the author and references sections.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dsuds
      integer iflag, il, ip, is, iwork(*), ks, kt, ku, kv, mlso, neq,
     1     nrda, nuk
      double precision a(nrda,*), b(*), work(*), x(*)
c
c***first executable statement  dsuds
      is = 2
      ip = 3
      il = ip + neq
      kv = 1 + neq
      kt = kv + neq
      ks = kt + neq
      ku = ks + nuk
c
      call dlssud(a,x,b,neq,nuk,nrda,work(ku),nuk,iflag,mlso,iwork(1),
     1            iwork(is),a,work(1),iwork(ip),b,work(kv),work(kt),
     2            iwork(il),work(ks))
c
      return
      end
