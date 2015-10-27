*deck dlssud
      subroutine dlssud (a, x, b, n, m, nrda, u, nrdu, iflag, mlso,
     +   irank, iscale, q, diag, kpivot, s, div, td, isflg, scales)
c***begin prologue  dlssud
c***subsidiary
c***purpose  subsidiary to dbvsup and dsuds
c***library   slatec
c***type      double precision (lssuds-s, dlssud-d)
c***author  watts, h. a., (snla)
c***description
c
c    dlssud solves the underdetermined system of equations  a z = b,
c    where a is n by m and n .le. m.  in particular, if rank a equals
c    ira, a vector x and a matrix u are determined such that x is the
c    unique solution of smallest length, satisfying a x = b, and the
c    columns of u form an orthonormal basis for the null space of a,
c    satisfying a u = 0 .  then all solutions z are given by
c              z = x + c(1)*u(1) + ..... + c(m-ira)*u(m-ira)
c    where u(j) represents the j-th column of u and the c(j) are
c    arbitrary constants.
c    if the system of equations are not compatible, only the least
c    squares solution of minimal length is computed.
c
c *********************************************************************
c   input
c *********************************************************************
c
c     a -- contains the matrix of n equations in m unknowns, a remains
c          unchanged, must be dimensioned nrda by m.
c     x -- solution array of length at least m.
c     b -- given constant vector of length n, b remains unchanged.
c     n -- number of equations, n greater or equal to 1.
c     m -- number of unknowns, m greater or equal to n.
c     nrda -- row dimension of a, nrda greater or equal to n.
c     u -- matrix used for solution, must be dimensioned nrdu by
c          (m - rank of a).
c          (storage for u may be ignored when only the minimal length
c           solution x is desired)
c     nrdu -- row dimension of u, nrdu greater or equal to m.
c             (if only the minimal length solution is wanted,
c              nrdu=0 is acceptable)
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
c                linear space defined by the matrix u.
c     irank -- variable used for the rank of a, set by the code.
c     iscale -- scaling indicator
c               =-1 if the matrix a is to be pre-scaled by
c               columns when appropriate.
c               if the scaling indicator is not equal to -1
c               no scaling will be attempted.
c            for most problems scaling will probably not be necessary.
c     q -- matrix used for the transformation, must be dimensioned
c            nrda by m.
c     diag,kpivot,s, -- arrays of length at least n used for internal
c      div,td,scales    storage (except for scales which is m).
c     isflg -- storage for an internal variable.
c
c *********************************************************************
c   output
c *********************************************************************
c
c     iflag -- status indicator
c            =1 if solution was obtained.
c            =2 if improper input is detected.
c            =3 if rank of matrix is less than n.
c               to continue, simply reset iflag=1 and call dlssud again.
c            =4 if the system of equations appears to be inconsistent.
c               however, the least squares solution of minimal length
c               was obtained.
c     x -- minimal length least squares solution of a z = b
c     irank -- numerically determined rank of a, must not be altered
c              on succeeding calls with input values of iflag=1.
c     u -- matrix whose m-irank columns are mutually orthogonal unit
c          vectors which span the null space of a. this is to be ignored
c          when mlso was set to zero or iflag=4 on output.
c     q -- contains the strictly upper triangular part of the reduced
c           matrix and transformation information.
c     diag -- contains the diagonal elements of the triangular reduced
c             matrix.
c     kpivot -- contains the pivotal information.  the row interchanges
c               performed on the original matrix are recorded here.
c     s -- contains the solution of the lower triangular system.
c     div,td -- contains transformation information for rank
c               deficient problems.
c     scales -- contains the column scaling parameters.
c
c *********************************************************************
c
c***see also  dbvsup, dsuds
c***references  h. a. watts, solving linear least squares problems
c                 using sods/suds/cods, sandia report sand77-0683,
c                 sandia laboratories, 1977.
c***routines called  d1mach, ddot, dohtrl, dorthr, j4save, xermax,
c                    xermsg, xgetf, xsetf
c***revision history  (yymmdd)
c   750601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900328  added type section.  (wrb)
c   910408  updated the author and references sections.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dlssud
      integer j4save
      double precision ddot, d1mach
      integer i, iflag, irank, irp, iscale, isflg, j, jr, k, kp,
     1     kpivot(*), l, m, maxmes, mj, mlso, n, nfat, nfatal, nmir,
     2     nrda, nrdu, nu
      double precision a(nrda,*), b(*), diag(*), div(*), gam, gamma,
     1     q(nrda,*), res, s(*), scales(*), ss, td(*), u(nrdu,*), uro,
     2     x(*)
c
c     ******************************************************************
c
c          machine precision (computer unit roundoff value) is defined
c          by the function d1mach.
c
c     ******************************************************************
c
c     begin block permitting ...exits to 310
c        begin block permitting ...exits to 80
c***first executable statement  dlssud
            uro = d1mach(4)
c
            if (n .lt. 1 .or. m .lt. n .or. nrda .lt. n) go to 70
            if (nrdu .ne. 0 .and. nrdu .lt. m) go to 70
               if (iflag .gt. 0) go to 60
c
                  call xgetf(nfatal)
                  maxmes = j4save(4,0,.false.)
                  isflg = -15
                  if (iflag .eq. 0) go to 10
                     isflg = iflag
                     nfat = -1
                     if (nfatal .eq. 0) nfat = 0
                     call xsetf(nfat)
                     call xermax(1)
   10             continue
c
c                 copy matrix a into matrix q
c
                  do 30 k = 1, m
                     do 20 j = 1, n
                        q(j,k) = a(j,k)
   20                continue
   30             continue
c
c                 use orthogonal transformations to reduce q to lower
c                 triangular form
c
                  call dorthr(q,n,m,nrda,iflag,irank,iscale,diag,kpivot,
     1                        scales,div,td)
c
                  call xsetf(nfatal)
                  call xermax(maxmes)
                  if (irank .eq. n) go to 40
c
c                    for rank deficient problems use additional
c                    orthogonal transformations to further reduce q
c
                     if (irank .ne. 0)
     1                  call dohtrl(q,n,nrda,diag,irank,div,td)
c     ...............exit
                     go to 310
   40             continue
c
c                 store divisors for the triangular solution
c
                  do 50 k = 1, n
                     div(k) = diag(k)
   50             continue
c        .........exit
                  go to 80
   60          continue
c        ......exit
               if (iflag .eq. 1) go to 80
   70       continue
c
c           invalid input for dlssud
            iflag = 2
            call xermsg ('slatec', 'dlssud',
     +         'invalid imput parameters.', 2, 1)
c     ......exit
            go to 310
   80    continue
c
c
         if (irank .gt. 0) go to 130
c
c           special case for the null matrix
            do 110 k = 1, m
               x(k) = 0.0d0
               if (mlso .eq. 0) go to 100
                  u(k,k) = 1.0d0
                  do 90 j = 1, m
                     if (j .ne. k) u(j,k) = 0.0d0
   90             continue
  100          continue
  110       continue
            do 120 k = 1, n
               if (b(k) .gt. 0.0d0) iflag = 4
  120       continue
         go to 300
  130    continue
c           begin block permitting ...exits to 180
c
c              copy constant vector into s after first interchanging
c              the elements according to the pivotal sequence
c
               do 140 k = 1, n
                  kp = kpivot(k)
                  x(k) = b(kp)
  140          continue
               do 150 k = 1, n
                  s(k) = x(k)
  150          continue
c
               irp = irank + 1
               nu = 1
               if (mlso .eq. 0) nu = 0
c           ...exit
               if (irank .eq. n) go to 180
c
c              for rank deficient problems we must apply the
c              orthogonal transformation to s
c              we also check to see if the system appears to be
c              inconsistent
c
               nmir = n - irank
               ss = ddot(n,s(1),1,s(1),1)
               do 170 l = 1, irank
                  k = irp - l
                  gam = ((td(k)*s(k)) + ddot(nmir,q(irp,k),1,s(irp),1))
     1                  /(td(k)*div(k))
                  s(k) = s(k) + gam*td(k)
                  do 160 j = irp, n
                     s(j) = s(j) + gam*q(j,k)
  160             continue
  170          continue
               res = ddot(nmir,s(irp),1,s(irp),1)
c           ...exit
               if (res
     1             .le. ss*(10.0d0*max(10.0d0**isflg,10.0d0*uro))**2)
     2            go to 180
c
c              inconsistent system
               iflag = 4
               nu = 0
  180       continue
c
c           apply forward substitution to solve lower triangular system
c
            s(1) = s(1)/div(1)
            if (irank .lt. 2) go to 200
            do 190 k = 2, irank
               s(k) = (s(k) - ddot(k-1,q(k,1),nrda,s(1),1))/div(k)
  190       continue
  200       continue
c
c           initialize x vector and then apply orthogonal transformation
c
            do 210 k = 1, m
               x(k) = 0.0d0
               if (k .le. irank) x(k) = s(k)
  210       continue
c
            do 230 jr = 1, irank
               j = irp - jr
               mj = m - j + 1
               gamma = ddot(mj,q(j,j),nrda,x(j),1)/(diag(j)*q(j,j))
               do 220 k = j, m
                  x(k) = x(k) + gamma*q(j,k)
  220          continue
  230       continue
c
c           rescale answers as dictated
c
            do 240 k = 1, m
               x(k) = x(k)*scales(k)
  240       continue
c
            if (nu .eq. 0 .or. m .eq. irank) go to 290
c
c              initialize u matrix and then apply orthogonal
c              transformation
c
               l = m - irank
               do 280 k = 1, l
                  do 250 i = 1, m
                     u(i,k) = 0.0d0
                     if (i .eq. irank + k) u(i,k) = 1.0d0
  250             continue
c
                  do 270 jr = 1, irank
                     j = irp - jr
                     mj = m - j + 1
                     gamma = ddot(mj,q(j,j),nrda,u(j,k),1)
     1                       /(diag(j)*q(j,j))
                     do 260 i = j, m
                        u(i,k) = u(i,k) + gamma*q(j,i)
  260                continue
  270             continue
  280          continue
  290       continue
  300    continue
  310 continue
c
      return
      end
