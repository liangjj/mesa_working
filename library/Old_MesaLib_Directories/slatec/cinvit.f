*deck cinvit
      subroutine cinvit (nm, n, ar, ai, wr, wi, select, mm, m, zr, zi,
     +   ierr, rm1, rm2, rv1, rv2)
c***begin prologue  cinvit
c***purpose  compute the eigenvectors of a complex upper hessenberg
c            associated with specified eigenvalues using inverse
c            iteration.
c***library   slatec (eispack)
c***category  d4c2b
c***type      complex (invit-s, cinvit-c)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of the algol procedure cxinvit
c     by peters and wilkinson.
c     handbook for auto. comp. vol.ii-linear algebra, 418-439(1971).
c
c     this subroutine finds those eigenvectors of a complex upper
c     hessenberg matrix corresponding to specified eigenvalues,
c     using inverse iteration.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameters, ar, ai, zr and zi, as declared in the
c          calling program dimension statement.  nm is an integer
c          variable.
c
c        n is the order of the matrix a=(ar,ai).  n is an integer
c          variable.  n must be less than or equal to nm.
c
c        ar and ai contain the real and imaginary parts, respectively,
c          of the complex upper hessenberg matrix.  ar and ai are
c          two-dimensional real arrays, dimensioned ar(nm,n)
c          and ai(nm,n).
c
c        wr and wi contain the real and imaginary parts, respectively,
c          of the eigenvalues of the matrix.  the eigenvalues must be
c          stored in a manner identical to that of subroutine  comlr,
c          which recognizes possible splitting of the matrix.  wr and
c          wi are one-dimensional real arrays, dimensioned wr(n) and
c          wi(n).
c
c        select specifies the eigenvectors to be found.  the
c          eigenvector corresponding to the j-th eigenvalue is
c          specified by setting select(j) to .true.  select is a
c          one-dimensional logical array, dimensioned select(n).
c
c        mm should be set to an upper bound for the number of
c          eigenvectors to be found.  mm is an integer variable.
c
c     on output
c
c        ar, ai, wi, and select are unaltered.
c
c        wr may have been altered since close eigenvalues are perturbed
c          slightly in searching for independent eigenvectors.
c
c        m is the number of eigenvectors actually found.  m is an
c          integer variable.
c
c        zr and zi contain the real and imaginary parts, respectively,
c          of the eigenvectors corresponding to the flagged eigenvalues.
c          the eigenvectors are normalized so that the component of
c          largest magnitude is 1.  any vector which fails the
c          acceptance test is set to zero.  zr and zi are
c          two-dimensional real arrays, dimensioned zr(nm,mm) and
c          zi(nm,mm).
c
c        ierr is an integer flag set to
c          zero       for normal return,
c          -(2*n+1)   if more than mm eigenvectors have been requested
c                     (the mm eigenvectors calculated to this point are
c                     in zr and zi),
c          -k         if the iteration corresponding to the k-th
c                     value fails (if this occurs more than once, k
c                     is the index of the last occurrence); the
c                     corresponding columns of zr and zi are set to
c                     zero vectors,
c          -(n+k)     if both error situations occur.
c
c        rv1 and rv2 are one-dimensional real arrays used for
c          temporary storage, dimensioned rv1(n) and rv2(n).
c          they hold the approximate eigenvectors during the inverse
c          iteration process.
c
c        rm1 and rm2 are two-dimensional real arrays used for
c          temporary storage, dimensioned rm1(n,n) and rm2(n,n).
c          these arrays hold the triangularized form of the upper
c          hessenberg matrix used in the inverse iteration process.
c
c     the algol procedure guessvec appears in cinvit in-line.
c
c     calls pythag(a,b) for sqrt(a**2 + b**2).
c     calls cdiv for complex division.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c     ------------------------------------------------------------------
c
c***references  b. t. smith, j. m. boyle, j. j. dongarra, b. s. garbow,
c                 y. ikebe, v. c. klema and c. b. moler, matrix eigen-
c                 system routines - eispack guide, springer-verlag,
c                 1976.
c***routines called  cdiv, pythag
c***revision history  (yymmdd)
c   760101  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  cinvit
c
      integer i,j,k,m,n,s,ii,mm,mp,nm,uk,ip1,its,km1,ierr
      real ar(nm,*),ai(nm,*),wr(*),wi(*),zr(nm,*),zi(nm,*)
      real rm1(n,*),rm2(n,*),rv1(*),rv2(*)
      real x,y,eps3,norm,normv,growto,ilambd,rlambd,ukroot
      real pythag
      logical select(n)
c
c***first executable statement  cinvit
      ierr = 0
      uk = 0
      s = 1
c
      do 980 k = 1, n
         if (.not. select(k)) go to 980
         if (s .gt. mm) go to 1000
         if (uk .ge. k) go to 200
c     .......... check for possible splitting ..........
         do 120 uk = k, n
            if (uk .eq. n) go to 140
            if (ar(uk+1,uk) .eq. 0.0e0 .and. ai(uk+1,uk) .eq. 0.0e0)
     1         go to 140
  120    continue
c     .......... compute infinity norm of leading uk by uk
c                (hessenberg) matrix ..........
  140    norm = 0.0e0
         mp = 1
c
         do 180 i = 1, uk
            x = 0.0e0
c
            do 160 j = mp, uk
  160       x = x + pythag(ar(i,j),ai(i,j))
c
            if (x .gt. norm) norm = x
            mp = i
  180    continue
c     .......... eps3 replaces zero pivot in decomposition
c                and close roots are modified by eps3 ..........
         if (norm .eq. 0.0e0) norm = 1.0e0
         eps3 = norm
  190    eps3 = 0.5e0*eps3
         if (norm + eps3 .gt. norm) go to 190
         eps3 = 2.0e0*eps3
c     .......... growto is the criterion for growth ..........
         ukroot = sqrt(real(uk))
         growto = 0.1e0 / ukroot
  200    rlambd = wr(k)
         ilambd = wi(k)
         if (k .eq. 1) go to 280
         km1 = k - 1
         go to 240
c     .......... perturb eigenvalue if it is close
c                to any previous eigenvalue ..........
  220    rlambd = rlambd + eps3
c     .......... for i=k-1 step -1 until 1 do -- ..........
  240    do 260 ii = 1, km1
            i = k - ii
            if (select(i) .and. abs(wr(i)-rlambd) .lt. eps3 .and.
     1         abs(wi(i)-ilambd) .lt. eps3) go to 220
  260    continue
c
         wr(k) = rlambd
c     .......... form upper hessenberg (ar,ai)-(rlambd,ilambd)*i
c                and initial complex vector ..........
  280    mp = 1
c
         do 320 i = 1, uk
c
            do 300 j = mp, uk
               rm1(i,j) = ar(i,j)
               rm2(i,j) = ai(i,j)
  300       continue
c
            rm1(i,i) = rm1(i,i) - rlambd
            rm2(i,i) = rm2(i,i) - ilambd
            mp = i
            rv1(i) = eps3
  320    continue
c     .......... triangular decomposition with interchanges,
c                replacing zero pivots by eps3 ..........
         if (uk .eq. 1) go to 420
c
         do 400 i = 2, uk
            mp = i - 1
            if (pythag(rm1(i,mp),rm2(i,mp)) .le.
     1         pythag(rm1(mp,mp),rm2(mp,mp))) go to 360
c
            do 340 j = mp, uk
               y = rm1(i,j)
               rm1(i,j) = rm1(mp,j)
               rm1(mp,j) = y
               y = rm2(i,j)
               rm2(i,j) = rm2(mp,j)
               rm2(mp,j) = y
  340       continue
c
  360       if (rm1(mp,mp) .eq. 0.0e0 .and. rm2(mp,mp) .eq. 0.0e0)
     1         rm1(mp,mp) = eps3
            call cdiv(rm1(i,mp),rm2(i,mp),rm1(mp,mp),rm2(mp,mp),x,y)
            if (x .eq. 0.0e0 .and. y .eq. 0.0e0) go to 400
c
            do 380 j = i, uk
               rm1(i,j) = rm1(i,j) - x * rm1(mp,j) + y * rm2(mp,j)
               rm2(i,j) = rm2(i,j) - x * rm2(mp,j) - y * rm1(mp,j)
  380       continue
c
  400    continue
c
  420    if (rm1(uk,uk) .eq. 0.0e0 .and. rm2(uk,uk) .eq. 0.0e0)
     1      rm1(uk,uk) = eps3
         its = 0
c     .......... back substitution
c                for i=uk step -1 until 1 do -- ..........
  660    do 720 ii = 1, uk
            i = uk + 1 - ii
            x = rv1(i)
            y = 0.0e0
            if (i .eq. uk) go to 700
            ip1 = i + 1
c
            do 680 j = ip1, uk
               x = x - rm1(i,j) * rv1(j) + rm2(i,j) * rv2(j)
               y = y - rm1(i,j) * rv2(j) - rm2(i,j) * rv1(j)
  680       continue
c
  700       call cdiv(x,y,rm1(i,i),rm2(i,i),rv1(i),rv2(i))
  720    continue
c     .......... acceptance test for eigenvector
c                and normalization ..........
         its = its + 1
         norm = 0.0e0
         normv = 0.0e0
c
         do 780 i = 1, uk
            x = pythag(rv1(i),rv2(i))
            if (normv .ge. x) go to 760
            normv = x
            j = i
  760       norm = norm + x
  780    continue
c
         if (norm .lt. growto) go to 840
c     .......... accept vector ..........
         x = rv1(j)
         y = rv2(j)
c
         do 820 i = 1, uk
            call cdiv(rv1(i),rv2(i),x,y,zr(i,s),zi(i,s))
  820    continue
c
         if (uk .eq. n) go to 940
         j = uk + 1
         go to 900
c     .......... in-line procedure for choosing
c                a new starting vector ..........
  840    if (its .ge. uk) go to 880
         x = ukroot
         y = eps3 / (x + 1.0e0)
         rv1(1) = eps3
c
         do 860 i = 2, uk
  860    rv1(i) = y
c
         j = uk - its + 1
         rv1(j) = rv1(j) - eps3 * x
         go to 660
c     .......... set error -- unaccepted eigenvector ..........
  880    j = 1
         ierr = -k
c     .......... set remaining vector components to zero ..........
  900    do 920 i = j, n
            zr(i,s) = 0.0e0
            zi(i,s) = 0.0e0
  920    continue
c
  940    s = s + 1
  980 continue
c
      go to 1001
c     .......... set error -- underestimate of eigenvector
c                space required ..........
 1000 if (ierr .ne. 0) ierr = ierr - n
      if (ierr .eq. 0) ierr = -(2 * n + 1)
 1001 m = s - 1
      return
      end
