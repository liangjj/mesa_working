*deck qrfac
      subroutine qrfac (m, n, a, lda, pivot, ipvt, lipvt, sigma, acnorm,
     +   wa)
c***begin prologue  qrfac
c***subsidiary
c***purpose  subsidiary to snls1, snls1e, snsq and snsqe
c***library   slatec
c***type      single precision (qrfac-s, dqrfac-d)
c***author  (unknown)
c***description
c
c     this subroutine uses householder transformations with column
c     pivoting (optional) to compute a qr factorization of the
c     m by n matrix a. that is, qrfac determines an orthogonal
c     matrix q, a permutation matrix p, and an upper trapezoidal
c     matrix r with diagonal elements of nonincreasing magnitude,
c     such that a*p = q*r. the householder transformation for
c     column k, k = 1,2,...,min(m,n), is of the form
c
c                           t
c           i - (1/u(k))*u*u
c
c     where u has zeros in the first k-1 positions. the form of
c     this transformation and the method of pivoting first
c     appeared in the corresponding linpack subroutine.
c
c     the subroutine statement is
c
c       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,sigma,acnorm,wa)
c
c     where
c
c       m is a positive integer input variable set to the number
c         of rows of a.
c
c       n is a positive integer input variable set to the number
c         of columns of a.
c
c       a is an m by n array. on input a contains the matrix for
c         which the qr factorization is to be computed. on output
c         the strict upper trapezoidal part of a contains the strict
c         upper trapezoidal part of r, and the lower trapezoidal
c         part of a contains a factored form of q (the non-trivial
c         elements of the u vectors described above).
c
c       lda is a positive integer input variable not less than m
c         which specifies the leading dimension of the array a.
c
c       pivot is a logical input variable. if pivot is set .true.,
c         then column pivoting is enforced. if pivot is set .false.,
c         then no column pivoting is done.
c
c       ipvt is an integer output array of length lipvt. ipvt
c         defines the permutation matrix p such that a*p = q*r.
c         column j of p is column ipvt(j) of the identity matrix.
c         if pivot is .false., ipvt is not referenced.
c
c       lipvt is a positive integer input variable. if pivot is
c             .false., then lipvt may be as small as 1. if pivot is
c             .true., then lipvt must be at least n.
c
c       sigma is an output array of length n which contains the
c         diagonal elements of r.
c
c       acnorm is an output array of length n which contains the
c         norms of the corresponding columns of the input matrix a.
c         if this information is not needed, then acnorm can coincide
c         with sigma.
c
c       wa is a work array of length n. if pivot is .false., then wa
c         can coincide with sigma.
c
c***see also  snls1, snls1e, snsq, snsqe
c***routines called  enorm, r1mach
c***revision history  (yymmdd)
c   800301  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900326  removed duplicate information from description section.
c           (wrb)
c   900328  added type section.  (wrb)
c***end prologue  qrfac
      integer m,n,lda,lipvt
      integer ipvt(*)
      logical pivot
      real a(lda,*),sigma(*),acnorm(*),wa(*)
      integer i,j,jp1,k,kmax,minmn
      real ajnorm,epsmch,one,p05,sum,temp,zero
      real r1mach,enorm
      save one, p05, zero
      data one,p05,zero /1.0e0,5.0e-2,0.0e0/
c***first executable statement  qrfac
      epsmch = r1mach(4)
c
c     compute the initial column norms and initialize several arrays.
c
      do 10 j = 1, n
         acnorm(j) = enorm(m,a(1,j))
         sigma(j) = acnorm(j)
         wa(j) = sigma(j)
         if (pivot) ipvt(j) = j
   10    continue
c
c     reduce a to r with householder transformations.
c
      minmn = min(m,n)
      do 110 j = 1, minmn
         if (.not.pivot) go to 40
c
c        bring the column of largest norm into the pivot position.
c
         kmax = j
         do 20 k = j, n
            if (sigma(k) .gt. sigma(kmax)) kmax = k
   20       continue
         if (kmax .eq. j) go to 40
         do 30 i = 1, m
            temp = a(i,j)
            a(i,j) = a(i,kmax)
            a(i,kmax) = temp
   30       continue
         sigma(kmax) = sigma(j)
         wa(kmax) = wa(j)
         k = ipvt(j)
         ipvt(j) = ipvt(kmax)
         ipvt(kmax) = k
   40    continue
c
c        compute the householder transformation to reduce the
c        j-th column of a to a multiple of the j-th unit vector.
c
         ajnorm = enorm(m-j+1,a(j,j))
         if (ajnorm .eq. zero) go to 100
         if (a(j,j) .lt. zero) ajnorm = -ajnorm
         do 50 i = j, m
            a(i,j) = a(i,j)/ajnorm
   50       continue
         a(j,j) = a(j,j) + one
c
c        apply the transformation to the remaining columns
c        and update the norms.
c
         jp1 = j + 1
         if (n .lt. jp1) go to 100
         do 90 k = jp1, n
            sum = zero
            do 60 i = j, m
               sum = sum + a(i,j)*a(i,k)
   60          continue
            temp = sum/a(j,j)
            do 70 i = j, m
               a(i,k) = a(i,k) - temp*a(i,j)
   70          continue
            if (.not.pivot .or. sigma(k) .eq. zero) go to 80
            temp = a(j,k)/sigma(k)
            sigma(k) = sigma(k)*sqrt(max(zero,one-temp**2))
            if (p05*(sigma(k)/wa(k))**2 .gt. epsmch) go to 80
            sigma(k) = enorm(m-j,a(jp1,k))
            wa(k) = sigma(k)
   80       continue
   90       continue
  100    continue
         sigma(j) = -ajnorm
  110    continue
      return
c
c     last card of subroutine qrfac.
c
      end
