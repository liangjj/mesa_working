*deck @(#)dsym.f	5.1  11/6/94
      subroutine dsym(triang,u,t2,t3,t4,nbf,nsym,nobs,
     #                csym,c)
      implicit integer(a-z)
c
c
      real*8 triang(*),u(*),t2(*),t3(nbf,nbf),t4(*)
      real*8 csym(nbf,nbf),c(nbf,nbf)
      real*8 big
      integer nobs(nsym)
c
c
      parameter (big=1.0d+44)
c
c ----  pack the fock operator
c
      kj=0
      jt=0
      jof=0
      do 1 i=1,nsym
         nob=nobs(i)
         do 2 j=1,nob
            do 3 k=1,j
               kj=kj+1
               u(kj)=triang(jt+k)
    3       continue
            jt=jt+jof+j
    2    continue
         jof=jof+nob
         jt=jt+nob
    1 continue
c
      call scopy(kj,u,1,triang,1)
c
      jx=1
      ixx=1
      ix=1
      do 11 i=1,nsym
         nob=nobs(i)
         nnp=nob*(nob+1)/2
         call degrsp(nob,nnp,triang(ix),t4(jx),1,u(ixx),t2,t3)
         call ebc(c(1,jx),csym(1,jx),u(ixx),nbf,nob,nob)
         ix=ix+nnp
         ixx=ixx+nob*nob
         jx=jx+nob
  11  continue
c
c ---- reorder the eigenvectors and eigenvalues into nosym order
c
      do 12 i=1,nbf
         jtest=ismin(nbf,t4,1)
         t2(i)=t4(jtest)
         t4(jtest)=big
         call scopy(nbf,c(1,jtest),1,t3(1,i),1)
  12  continue
c
      call scopy(nbf,t2,1,t4,1)
      call scopy(nbf*nbf,t3,1,c,1)
c

      return
      end
