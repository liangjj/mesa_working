*deck dmgsbv
      subroutine dmgsbv (m, n, a, ia, niv, iflag, s, p, ip, inhomo, v,
     +   w, wcnd)
c***begin prologue  dmgsbv
c***subsidiary
c***purpose  subsidiary to dbvsup
c***library   slatec
c***type      double precision (mgsbv-s, dmgsbv-d)
c***author  watts, h. a., (snla)
c***description
c
c **********************************************************************
c orthogonalize a set of n double precision vectors and determine their
c rank.
c
c **********************************************************************
c input
c **********************************************************************
c   m = dimension of vectors.
c   n = no. of vectors.
c   a = array whose first n cols contain the vectors.
c   ia = first dimension of array a (col length).
c   niv = number of independent vectors needed.
c   inhomo = 1 corresponds to having a non-zero particular solution.
c   v = particular solution vector (not included in the pivoting).
c   indpvt = 1 means pivoting will not be used.
c
c **********************************************************************
c output
c **********************************************************************
c   niv = no. of linear independent vectors in input set.
c     a = matrix whose first niv cols. contain niv orthogonal vectors
c         which span the vector space determined by the input vectors.
c   iflag
c          = 0 success
c          = 1 incorrect input
c          = 2 rank of new vectors less than n
c   p = decomposition matrix.  p is upper triangular and
c             (old vectors) = (new vectors) * p.
c         the old vectors will be reordered due to pivoting.
c         the dimension of p must be .ge. n*(n+1)/2.
c             (  n*(2*n+1) when n .ne. nfcc )
c   ip = pivoting vector. the dimension of ip must be .ge. n.
c             (  2*n when n .ne. nfcc )
c   s = square of norms of incoming vectors.
c   v = vector which is orthogonal to the vectors of a.
c   w = orthogonalization information for the vector v.
c   wcnd = worst case (smallest) norm decrement value of the
c          vectors being orthogonalized  (represents a test
c          for linear dependence of the vectors).
c **********************************************************************
c
c***see also  dbvsup
c***routines called  ddot, dprvec
c***common blocks    dml18j, dml5mc
c***revision history  (yymmdd)
c   750601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890921  realigned order of variables in certain common blocks.
c           (wrb)
c   890921  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  dmgsbv
c
      double precision ddot, dprvec
      integer i, ia, icoco, iflag, indpvt, inhomo, integ, ip(*), ip1,
     1     ix, iz, j, jk, jp, jq, jy, jz, k, kd, kj, kp, l, lix, lpar,
     2     lr, m, m2, mxnon, n, ndisk, neq, neqivp, nfcc, nic, niv,
     3     nivn, nmnr, nn, nopg, np1, nps, nr, nrm1, ntape, ntp,
     4     numort, nxpts
      double precision a(ia,*), ae, dot, eps, fouru, p(*), pjp, psave,
     1     re, ry, s(*), sqovfl, sru, sv, t, tol, twou, uro, v(*), vl,
     2     vnorm, w(*), wcnd, y
c
c
      common /dml18j/ ae,re,tol,nxpts,nic,nopg,mxnon,ndisk,ntape,neq,
     1                indpvt,integ,nps,ntp,neqivp,numort,nfcc,
     2                icoco
c
      common /dml5mc/ uro,sru,eps,sqovfl,twou,fouru,lpar
c
c***first executable statement  dmgsbv
      if (m .gt. 0 .and. n .gt. 0 .and. ia .ge. m) go to 10
         iflag = 1
      go to 280
   10 continue
c        begin block permitting ...exits to 270
c           begin block permitting ...exits to 260
c
               jp = 0
               iflag = 0
               np1 = n + 1
               y = 0.0d0
               m2 = m/2
c
c              calculate square of norms of incoming vectors and search
c              for vector with largest magnitude
c
               j = 0
               do 40 i = 1, n
                  vl = ddot(m,a(1,i),1,a(1,i),1)
                  s(i) = vl
                  if (n .eq. nfcc) go to 20
                     j = 2*i - 1
                     p(j) = vl
                     ip(j) = j
   20             continue
                  j = j + 1
                  p(j) = vl
                  ip(j) = j
                  if (vl .le. y) go to 30
                     y = vl
                     ix = i
   30             continue
   40          continue
               if (indpvt .ne. 1) go to 50
                  ix = 1
                  y = p(1)
   50          continue
               lix = ix
               if (n .ne. nfcc) lix = 2*ix - 1
               p(lix) = p(1)
               s(np1) = 0.0d0
               if (inhomo .eq. 1) s(np1) = ddot(m,v,1,v,1)
               wcnd = 1.0d0
               nivn = niv
               niv = 0
c
c           ...exit
               if (y .eq. 0.0d0) go to 260
c              *********************************************************
               do 240 nr = 1, n
c                 begin block permitting ...exits to 230
c              ......exit
                     if (nivn .eq. niv) go to 250
                     niv = nr
                     if (ix .eq. nr) go to 130
c
c                       pivoting of columns of p matrix
c
                        nn = n
                        lix = ix
                        lr = nr
                        if (n .eq. nfcc) go to 60
                           nn = nfcc
                           lix = 2*ix - 1
                           lr = 2*nr - 1
   60                   continue
                        if (nr .eq. 1) go to 80
                           kd = lix - lr
                           kj = lr
                           nrm1 = lr - 1
                           do 70 j = 1, nrm1
                              psave = p(kj)
                              jk = kj + kd
                              p(kj) = p(jk)
                              p(jk) = psave
                              kj = kj + nn - j
   70                      continue
                           jy = jk + nmnr
                           jz = jy - kd
                           p(jy) = p(jz)
   80                   continue
                        iz = ip(lix)
                        ip(lix) = ip(lr)
                        ip(lr) = iz
                        sv = s(ix)
                        s(ix) = s(nr)
                        s(nr) = sv
                        if (n .eq. nfcc) go to 110
                           if (nr .eq. 1) go to 100
                              kj = lr + 1
                              do 90 k = 1, nrm1
                                 psave = p(kj)
                                 jk = kj + kd
                                 p(kj) = p(jk)
                                 p(jk) = psave
                                 kj = kj + nfcc - k
   90                         continue
  100                      continue
                           iz = ip(lix+1)
                           ip(lix+1) = ip(lr+1)
                           ip(lr+1) = iz
  110                   continue
c
c                       pivoting of columns of vectors
c
                        do 120 l = 1, m
                           t = a(l,ix)
                           a(l,ix) = a(l,nr)
                           a(l,nr) = t
  120                   continue
  130                continue
c
c                    calculate p(nr,nr) as norm squared of pivotal
c                    vector
c
                     jp = jp + 1
                     p(jp) = y
                     ry = 1.0d0/y
                     nmnr = n - nr
                     if (n .eq. nfcc) go to 140
                        nmnr = nfcc - (2*nr - 1)
                        jp = jp + 1
                        p(jp) = 0.0d0
                        kp = jp + nmnr
                        p(kp) = y
  140                continue
                     if (nr .eq. n .or. nivn .eq. niv) go to 200
c
c                       calculate orthogonal projection vectors and
c                       search for largest norm
c
                        y = 0.0d0
                        ip1 = nr + 1
                        ix = ip1
c                       ************************************************
                        do 190 j = ip1, n
                           dot = ddot(m,a(1,nr),1,a(1,j),1)
                           jp = jp + 1
                           jq = jp + nmnr
                           if (n .ne. nfcc) jq = jq + nmnr - 1
                           p(jq) = p(jp) - dot*(dot*ry)
                           p(jp) = dot*ry
                           do 150 i = 1, m
                              a(i,j) = a(i,j) - p(jp)*a(i,nr)
  150                      continue
                           if (n .eq. nfcc) go to 170
                              kp = jp + nmnr
                              jp = jp + 1
                              pjp = ry*dprvec(m,a(1,nr),a(1,j))
                              p(jp) = pjp
                              p(kp) = -pjp
                              kp = kp + 1
                              p(kp) = ry*dot
                              do 160 k = 1, m2
                                 l = m2 + k
                                 a(k,j) = a(k,j) - pjp*a(l,nr)
                                 a(l,j) = a(l,j) + pjp*a(k,nr)
  160                         continue
                              p(jq) = p(jq) - pjp*(pjp/ry)
  170                      continue
c
c                          test for cancellation in recurrence relation
c
                           if (p(jq) .le. s(j)*sru)
     1                        p(jq) = ddot(m,a(1,j),1,a(1,j),1)
                           if (p(jq) .le. y) go to 180
                              y = p(jq)
                              ix = j
  180                      continue
  190                   continue
                        if (n .ne. nfcc) jp = kp
c                       ************************************************
                        if (indpvt .eq. 1) ix = ip1
c
c                       recompute norm squared of pivotal vector with
c                       scalar product
c
                        y = ddot(m,a(1,ix),1,a(1,ix),1)
c           ............exit
                        if (y .le. eps*s(ix)) go to 260
                        wcnd = min(wcnd,y/s(ix))
  200                continue
c
c                    compute orthogonal projection of particular
c                    solution
c
c                 ...exit
                     if (inhomo .ne. 1) go to 230
                     lr = nr
                     if (n .ne. nfcc) lr = 2*nr - 1
                     w(lr) = ddot(m,a(1,nr),1,v,1)*ry
                     do 210 i = 1, m
                        v(i) = v(i) - w(lr)*a(i,nr)
  210                continue
c                 ...exit
                     if (n .eq. nfcc) go to 230
                     lr = 2*nr
                     w(lr) = ry*dprvec(m,v,a(1,nr))
                     do 220 k = 1, m2
                        l = m2 + k
                        v(k) = v(k) + w(lr)*a(l,nr)
                        v(l) = v(l) - w(lr)*a(k,nr)
  220                continue
  230             continue
  240          continue
  250          continue
c              *********************************************************
c
c                  test for linear dependence of particular solution
c
c        ......exit
               if (inhomo .ne. 1) go to 270
               if ((n .gt. 1) .and. (s(np1) .lt. 1.0)) go to 270
               vnorm = ddot(m,v,1,v,1)
               if (s(np1) .ne. 0.0d0) wcnd = min(wcnd,vnorm/s(np1))
c        ......exit
               if (vnorm .ge. eps*s(np1)) go to 270
  260       continue
            iflag = 2
            wcnd = eps
  270    continue
  280 continue
      return
      end
