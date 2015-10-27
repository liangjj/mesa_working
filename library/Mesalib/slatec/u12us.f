*deck u12us
      subroutine u12us (a, mda, m, n, b, mdb, nb, mode, krank, rnorm, h,
     +   w, ir, ic)
c***begin prologue  u12us
c***subsidiary
c***purpose  subsidiary to ulsia
c***library   slatec
c***type      single precision (u12us-s, du12us-d)
c***author  (unknown)
c***description
c
c        given the householder lq factorization of a, this
c        subroutine solves the system ax=b. if the system
c        is of reduced rank, this routine returns a solution
c        according to the selected mode.
c
c       note - if mode.ne.2, w is never accessed.
c
c***see also  ulsia
c***routines called  saxpy, sdot, snrm2, sswap
c***revision history  (yymmdd)
c   810801  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c***end prologue  u12us
      dimension a(mda,*),b(mdb,*),rnorm(*),h(*),w(*)
      integer ic(*),ir(*)
c***first executable statement  u12us
      k=krank
      kp1=k+1
c
c        rank=0
c
      if(k.gt.0) go to 410
      do 404 jb=1,nb
      rnorm(jb)=snrm2(m,b(1,jb),1)
  404 continue
      do 406 jb=1,nb
      do 406 i=1,n
      b(i,jb)=0.0
  406 continue
      return
c
c     reorder b to reflect row interchanges
c
  410 continue
      i=0
  412 i=i+1
      if(i.eq.m) go to 418
      j=ir(i)
      if(j.eq.i) go to 412
      if(j.lt.0) go to 412
      ir(i)=-ir(i)
      do 413 jb=1,nb
      rnorm(jb)=b(i,jb)
  413 continue
      ij=i
  414 do 415 jb=1,nb
      b(ij,jb)=b(j,jb)
  415 continue
      ij=j
      j=ir(ij)
      ir(ij)=-ir(ij)
      if(j.ne.i) go to 414
      do 416 jb=1,nb
      b(ij,jb)=rnorm(jb)
  416 continue
      go to 412
  418 continue
      do 420 i=1,m
      ir(i)=abs(ir(i))
  420 continue
c
c     if a is of reduced rank and mode=2,
c     apply householder transformations to b
c
      if(mode.lt.2 .or. k.eq.m) go to 440
      mmk=m-k
      do 430 jb=1,nb
      do 425 j=1,k
      i=kp1-j
      tt=-sdot(mmk,a(kp1,i),1,b(kp1,jb),1)/w(i)
      tt=tt-b(i,jb)
      call saxpy(mmk,tt,a(kp1,i),1,b(kp1,jb),1)
      b(i,jb)=b(i,jb)+tt*w(i)
  425 continue
  430 continue
c
c     find norms of residual vector(s)..(before overwrite b)
c
  440 do 442 jb=1,nb
      rnorm(jb)=snrm2((m-k),b(kp1,jb),1)
  442 continue
c
c     back solve lower triangular l
c
      do 450 jb=1,nb
      do 448 i=1,k
      b(i,jb)=b(i,jb)/a(i,i)
      if(i.eq.k) go to 450
      ip1=i+1
      call saxpy(k-i,-b(i,jb),a(ip1,i),1,b(ip1,jb),1)
  448 continue
  450 continue
c
c
c      truncated solution
c
      if(k.eq.n) go to 462
      do 460 jb=1,nb
      do 460 i=kp1,n
      b(i,jb)=0.0
  460 continue
c
c     apply householder transformations to b
c
  462 do 470 i=1,k
      j=kp1-i
      tt=a(j,j)
      a(j,j)=h(j)
      do 465 jb=1,nb
      bb=-sdot(n-j+1,a(j,j),mda,b(j,jb),1)/h(j)
      call saxpy(n-j+1,bb,a(j,j),mda,b(j,jb),1)
  465 continue
      a(j,j)=tt
  470 continue
c
c
c     reorder b to reflect column interchanges
c
      i=0
  482 i=i+1
      if(i.eq.n) go to 488
      j=ic(i)
      if(j.eq.i) go to 482
      if(j.lt.0) go to 482
      ic(i)=-ic(i)
  484 call sswap(nb,b(j,1),mdb,b(i,1),mdb)
      ij=ic(j)
      ic(j)=-ic(j)
      j=ij
      if(j.eq.i) go to 482
      go to 484
  488 continue
      do 490 i=1,n
      ic(i)=abs(ic(i))
  490 continue
c
c        solution vectors are in first n rows of b(,)
c
      return
      end
