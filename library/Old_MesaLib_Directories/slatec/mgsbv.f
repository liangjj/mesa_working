*deck mgsbv
      subroutine mgsbv (m, n, a, ia, niv, iflag, s, p, ip, inhomo, v, w,
     +   wcnd)
c***begin prologue  mgsbv
c***subsidiary
c***purpose  subsidiary to bvsup
c***library   slatec
c***type      single precision (mgsbv-s, dmgsbv-d)
c***author  watts, h. a., (snla)
c***description
c
c **********************************************************************
c orthogonalize a set of n real vectors and determine their rank
c
c **********************************************************************
c input
c **********************************************************************
c   m = dimension of vectors
c   n = no. of vectors
c   a = array whose first n cols contain the vectors
c   ia = first dimension of array a (col length)
c   niv = number of independent vectors needed
c   inhomo = 1 corresponds to having a non-zero particular solution
c   v = particular solution vector (not included in the pivoting)
c   indpvt = 1 means pivoting will not be used
c
c **********************************************************************
c output
c **********************************************************************
c   niv = no. of linear independent vectors in input set
c     a = matrix whose first niv cols. contain niv orthogonal vectors
c         which span the vector space determined by the input vectors
c   iflag
c          = 0 success
c          = 1 incorrect input
c          = 2 rank of new vectors less than n
c   p = decomposition matrix.  p is upper triangular and
c             (old vectors) = (new vectors) * p.
c         the old vectors will be reordered due to pivoting
c         the dimension of p must be .ge. n*(n+1)/2.
c             (  n*(2*n+1) when n .ne. nfcc )
c   ip = pivoting vector. the dimension of ip must be .ge. n.
c             (  2*n when n .ne. nfcc )
c   s = square of norms of incoming vectors
c   v = vector which is orthogonal to the vectors of a
c   w = orthogonalization information for the vector v
c   wcnd = worst case (smallest) norm decrement value of the
c          vectors being orthogonalized  (represents a test
c          for linear dependence of the vectors)
c **********************************************************************
c
c***see also  bvsup
c***routines called  prvec, sdot
c***common blocks    ml18jr, ml5mco
c***revision history  (yymmdd)
c   750601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890921  realigned order of variables in certain common blocks.
c           (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  mgsbv
c
      dimension a(ia,*),v(*),w(*),p(*),ip(*),s(*)
c
c
      common /ml18jr/ ae,re,tol,nxpts,nic,nopg,mxnon,ndisk,ntape,neq,
     1                indpvt,integ,nps,ntp,neqivp,numort,nfcc,
     2                icoco
c
      common /ml5mco/ uro,sru,eps,sqovfl,twou,fouru,lpar
c
c***first executable statement  mgsbv
      if(m .gt. 0  .and.  n .gt. 0  .and.  ia .ge. m) go to 10
      iflag=1
      return
c
   10 jp=0
      iflag=0
      np1=n+1
      y=0.0
      m2=m/2
c
c     calculate square of norms of incoming vectors and search for
c     vector with largest magnitude
c
      j=0
      do 30 i=1,n
      vl=sdot(m,a(1,i),1,a(1,i),1)
      s(i)=vl
      if (n .eq. nfcc) go to 25
      j=2*i-1
      p(j)=vl
      ip(j)=j
   25 j=j+1
      p(j)=vl
      ip(j)=j
      if(vl .le. y) go to 30
      y=vl
      ix=i
   30 continue
      if (indpvt .ne. 1) go to 33
      ix=1
      y=p(1)
   33 lix=ix
      if (n .ne. nfcc) lix=2*ix-1
      p(lix)=p(1)
      s(np1)=0.
      if (inhomo .eq. 1) s(np1)=sdot(m,v,1,v,1)
      wcnd=1.
      nivn=niv
      niv=0
c
      if(y .eq. 0.0) go to 170
c **********************************************************************
      do 140 nr=1,n
      if (nivn .eq. niv) go to 150
      niv=nr
      if(ix .eq. nr) go to 80
c
c     pivoting of columns of p matrix
c
      nn=n
      lix=ix
      lr=nr
      if (n .eq. nfcc) go to 40
      nn=nfcc
      lix=2*ix-1
      lr=2*nr-1
   40 if(nr .eq. 1) go to 60
      kd=lix-lr
      kj=lr
      nrm1=lr-1
      do 50 j=1,nrm1
      psave=p(kj)
      jk=kj+kd
      p(kj)=p(jk)
      p(jk)=psave
   50 kj=kj+nn-j
      jy=jk+nmnr
      jz=jy-kd
      p(jy)=p(jz)
   60 iz=ip(lix)
      ip(lix)=ip(lr)
      ip(lr)=iz
      sv=s(ix)
      s(ix)=s(nr)
      s(nr)=sv
      if (n .eq. nfcc) go to 69
      if (nr .eq. 1) go to 67
      kj=lr+1
      do 65 k=1,nrm1
      psave=p(kj)
      jk=kj+kd
      p(kj)=p(jk)
      p(jk)=psave
   65 kj=kj+nfcc-k
   67 iz=ip(lix+1)
      ip(lix+1)=ip(lr+1)
      ip(lr+1)=iz
c
c     pivoting of columns of vectors
c
   69 do 70 l=1,m
      t=a(l,ix)
      a(l,ix)=a(l,nr)
   70 a(l,nr)=t
c
c     calculate p(nr,nr) as norm squared of pivotal vector
c
   80 jp=jp+1
      p(jp)=y
      ry=1.0/y
      nmnr=n-nr
      if (n .eq. nfcc) go to 85
      nmnr=nfcc-(2*nr-1)
      jp=jp+1
      p(jp)=0.
      kp=jp+nmnr
      p(kp)=y
   85 if(nr .eq. n  .or.  nivn .eq. niv) go to 125
c
c    calculate orthogonal projection vectors and search for largest norm
c
      y=0.0
      ip1=nr+1
      ix=ip1
c     ****************************************
      do 120 j=ip1,n
      dot=sdot(m,a(1,nr),1,a(1,j),1)
      jp=jp+1
      jq=jp+nmnr
      if (n .ne. nfcc) jq=jq+nmnr-1
      p(jq)=p(jp)-dot*(dot*ry)
      p(jp)=dot*ry
      do 90 i = 1,m
   90 a(i,j)=a(i,j)-p(jp)*a(i,nr)
      if (n .eq. nfcc) go to 99
      kp=jp+nmnr
      jp=jp+1
      pjp=ry*prvec(m,a(1,nr),a(1,j))
      p(jp)=pjp
      p(kp)=-pjp
      kp=kp+1
      p(kp)=ry*dot
      do 95 k=1,m2
      l=m2+k
      a(k,j)=a(k,j)-pjp*a(l,nr)
   95 a(l,j)=a(l,j)+pjp*a(k,nr)
      p(jq)=p(jq)-pjp*(pjp/ry)
c
c     test for cancellation in recurrence relation
c
   99 if(p(jq) .gt. s(j)*sru) go to 100
      p(jq)=sdot(m,a(1,j),1,a(1,j),1)
  100 if(p(jq) .le. y) go to 120
      y=p(jq)
      ix=j
  120 continue
      if (n .ne. nfcc) jp=kp
c     ****************************************
      if(indpvt .eq. 1) ix=ip1
c
c     recompute norm squared of pivotal vector with scalar product
c
      y=sdot(m,a(1,ix),1,a(1,ix),1)
      if(y  .le.  eps*s(ix))  go to 170
      wcnd=min(wcnd,y/s(ix))
c
c     compute orthogonal projection of particular solution
c
  125 if(inhomo .ne. 1) go to 140
      lr=nr
      if (n .ne. nfcc) lr=2*nr-1
      w(lr)=sdot(m,a(1,nr),1,v,1)*ry
      do 130 i=1,m
  130 v(i)=v(i)-w(lr)*a(i,nr)
      if (n .eq. nfcc) go to 140
      lr=2*nr
      w(lr)=ry*prvec(m,v,a(1,nr))
      do 135 k=1,m2
      l=m2+k
      v(k)=v(k)+w(lr)*a(l,nr)
  135 v(l)=v(l)-w(lr)*a(k,nr)
  140 continue
c **********************************************************************
c
c     test for linear dependence of particular solution
c
  150 if(inhomo .ne. 1) return
      if ((n .gt. 1) .and. (s(np1) .lt. 1.0)) return
      vnorm=sdot(m,v,1,v,1)
      if (s(np1) .ne. 0.) wcnd=min(wcnd,vnorm/s(np1))
      if(vnorm .ge. eps*s(np1)) return
  170 iflag=2
      wcnd=eps
      return
      end
