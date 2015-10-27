*deck cpqr79
      subroutine cpqr79 (ndeg, coeff, root, ierr, work)
c***begin prologue  cpqr79
c***purpose  find the zeros of a polynomial with complex coefficients.
c***library   slatec
c***category  f1a1b
c***type      complex (rpqr79-s, cpqr79-c)
c***keywords  complex polynomial, polynomial roots, polynomial zeros
c***author  vandevender, w. h., (snla)
c***description
c
c   abstract
c       this routine computes all zeros of a polynomial of degree ndeg
c       with complex coefficients by computing the eigenvalues of the
c       companion matrix.
c
c   description of parameters
c       the user must dimension all arrays appearing in the call list
c            coeff(ndeg+1), root(ndeg), work(2*ndeg*(ndeg+1))
c
c    --input--
c      ndeg    degree of polynomial
c
c      coeff   complex coefficients in descending order.  i.e.,
c              p(z)= coeff(1)*(z**ndeg) + coeff(ndeg)*z + coeff(ndeg+1)
c
c      work    real work array of dimension at least 2*ndeg*(ndeg+1)
c
c   --output--
c      root    complex vector of roots
c
c      ierr    output error code
c           - normal code
c          0  means the roots were computed.
c           - abnormal codes
c          1  more than 30 qr iterations on some eigenvalue of the
c             companion matrix
c          2  coeff(1)=0.0
c          3  ndeg is invalid (less than or equal to 0)
c
c***references  (none)
c***routines called  comqr, xermsg
c***revision history  (yymmdd)
c   791201  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c   911010  code reworked and simplified.  (rwc and wrb)
c***end prologue  cpqr79
      complex coeff(*), root(*), scale, c
      real work(*)
      integer ndeg, ierr, k, khr, khi, kwr, kwi, kad, kj
c***first executable statement  cpqr79
      ierr = 0
      if (abs(coeff(1)) .eq. 0.0) then
         ierr = 2
         call xermsg ('slatec', 'cpqr79',
     +      'leading coefficient is zero.', 2, 1)
         return
      endif
c
      if (ndeg .le. 0) then
         ierr = 3
         call xermsg ('slatec', 'cpqr79', 'degree invalid.', 3, 1)
         return
      endif
c
      if (ndeg .eq. 1) then
         root(1) = -coeff(2)/coeff(1)
         return
      endif
c
      scale = 1.0e0/coeff(1)
      khr = 1
      khi = khr+ndeg*ndeg
      kwr = khi+khi-khr
      kwi = kwr+ndeg
c
      do 10 k=1,kwr
         work(k) = 0.0e0
   10 continue
c
      do 20 k=1,ndeg
         kad = (k-1)*ndeg+1
         c = scale*coeff(k+1)
         work(kad) = -real(c)
         kj = khi+kad-1
         work(kj) = -aimag(c)
         if (k .ne. ndeg) work(kad+k) = 1.0e0
   20 continue
c
      call comqr (ndeg,ndeg,1,ndeg,work(khr),work(khi),work(kwr),
     1   work(kwi),ierr)
c
      if (ierr .ne. 0) then
         ierr = 1
         call xermsg ('slatec', 'cpqr79',
     +      'no convergence in 30 qr iterations.', 1, 1)
         return
      endif
c
      do 30 k=1,ndeg
         km1 = k-1
         root(k) = cmplx(work(kwr+km1),work(kwi+km1))
   30 continue
      return
      end
