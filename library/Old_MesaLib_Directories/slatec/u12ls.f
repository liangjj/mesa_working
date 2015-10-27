*deck u12ls
      subroutine u12ls (a, mda, m, n, b, mdb, nb, mode, krank, rnorm, h,
     +   w, ic, ir)
c***begin prologue  u12ls
c***subsidiary
c***purpose  subsidiary to llsia
c***library   slatec
c***type      single precision (u12ls-s, du12ls-d)
c***author  (unknown)
c***description
c
c        given the householder qr factorization of a, this
c        subroutine solves the system ax=b. if the system
c        is of reduced rank, this routine returns a solution
c        according to the selected mode.
c
c       note - if mode.ne.2, w is never accessed.
c
c***see also  llsia
c***routines called  saxpy, sdot, snrm2, sswap
c***revision history  (yymmdd)
c   810801  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c***end prologue  u12ls
      dimension a(mda,*),b(mdb,*),rnorm(*),h(*),w(*)
      integer ic(*),ir(*)
c***first executable statement  u12ls
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
c     apply householder transformations to b
c
      do 430 j=1,k
      tt=a(j,j)
      a(j,j)=h(j)
      do 425 i=1,nb
      bb=-sdot(m-j+1,a(j,j),1,b(j,i),1)/h(j)
      call saxpy(m-j+1,bb,a(j,j),1,b(j,i),1)
  425 continue
      a(j,j)=tt
  430 continue
c
c        find norms of residual vector(s)..(before overwrite b)
c
      do 440 jb=1,nb
      rnorm(jb)=snrm2((m-k),b(kp1,jb),1)
  440 continue
c
c     back solve upper triangular r
c
      i=k
  442 do 444 jb=1,nb
      b(i,jb)=b(i,jb)/a(i,i)
  444 continue
      if(i.eq.1) go to 450
      im1=i-1
      do 448 jb=1,nb
      call saxpy(im1,-b(i,jb),a(1,i),1,b(1,jb),1)
  448 continue
      i=im1
      go to 442
  450 continue
c
c     rank lt n
c
c      truncated solution
c
      if(k.eq.n) go to 480
      do 460 jb=1,nb
      do 460 i=kp1,n
      b(i,jb)=0.0
  460 continue
      if(mode.eq.1) go to 480
c
c      minimal length solution
c
      nmk=n-k
      do 470 jb=1,nb
      do 465 i=1,k
      tt=-sdot(nmk,a(i,kp1),mda,b(kp1,jb),1)/w(i)
      tt=tt-b(i,jb)
      call saxpy(nmk,tt,a(i,kp1),mda,b(kp1,jb),1)
      b(i,jb)=b(i,jb)+tt*w(i)
  465 continue
  470 continue
c
c
c     reorder b to reflect column interchanges
c
  480 continue
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
