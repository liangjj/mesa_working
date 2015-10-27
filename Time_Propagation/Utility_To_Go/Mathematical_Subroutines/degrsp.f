*deck @(#)degrsp.f	5.1  11/6/94
      subroutine degrsp(n,nnp,z,lambda,job,u,t1,t2)
c
c***begin prologue     degrsp.f
c***date written       860813  (yymmdd)
c***revision date      11/6/94
c   august 3, 1994     rlm at lanl
c      modifying the definition of the reldif statement so that it
c      performs in the unusual event of degenerate zero eigenvalues.
c***keywords           matrix, diagonalization, degenerate eigenvalues
c***author             saxe, paul (lanl)
c***source
c***purpose            diagonalization of a symmetric-packed matrix with
c                      rotation of degenerate eigenvectors to sensible
c                      arrangement.
c***description
c                      call degrsp(n,nnp,z,lambda,job,u,t1,t2)
c
c      n         integer
c                the matrix dimension
c
c      nnp       integer
c                the size of the packed matrix [n*(n+1)/2]
c
c      z         real (nnp)
c                the packed symmetric matrix to diagonalize
c
c      lambda    real (n)
c                the eigenvalues returned.
c
c      job       integer
c                0 if only eigenvalues needed,
c                1 if eigenvectors also needed
c
c      u         real (n,n)
c                the array for eigenvectors if job=1
c
c      t1, t2    real (n)
c                scratch arrays.
c
c      this routine currently handles doubly and triply degenerate eigenvalues
c
c***references
c***routines called    rsp(clams)
c***end prologue       degrsp
c
      implicit integer (a-z)
c
      real*8 z(nnp),u(n,n),lambda(n),t1(n),t2(n)
      real*8 fuzz,reldif,dif,a,b,small,sdot
      common/io/inp,iout
c
      parameter (fuzz=1.0d-08,small=1.0d-01)
c
c     ----- statement functions -----
c
      reldif(a,b)=abs((a-b)/(a+b))
      dif(a,b)=abs(a-b)
c
c     ----- diagonalize a and get eigenvalues and vectors -----
c
      call rsp(n,n,nnp,z,lambda,job,u,t1,t2,ierr)
c
      if (ierr.ne.0) call lnkerr('rsp returned non-zero error')
c
c     ----- if only eigenvalues needed, return -----
c
      if (job.eq.0) return
c
c----------------------------------------------------------------------
c   rotate eigenvectors so that they are unique. strategy is to
c   subtract the first of a degenerate pair from the second with a
c   coefficient such that the first largish element of the second
c   eigenvector is annihilated. the second vector is then projected from
c   the first to hopefully leave two nice vectors.
c-----------------------------------------------------------------------
c
c     ----- handle doubly degenerate roots for now -----
c
      i=1
    1 continue
c
c        --- replacing the relative difference with a direct absolute
c            difference.  this should still spot degeneracies.
c        if (reldif(lambda(i),lambda(i+1)).lt.fuzz) then
         if (dif(lambda(i),lambda(i+1)).lt.fuzz) then
            if (i.lt.n-1) then
c              if (reldif(lambda(i+1),lambda(i+2)).lt.fuzz) then
               if (dif(lambda(i+1),lambda(i+2)).lt.fuzz) then
c
c                 ----- first part for triply degenerate roots -----
c                       find reasonably large numbers to work with
c
                  do 2 j=1,n
                     if (abs(u(j,i))+abs(u(j,i+1)).gt.small) go to 3
    2             continue
                  call lnkerr('no value larger than small ????')
    3             continue
c
c                 ----- reorder so largest value is first -----
c
                  if (abs(u(j,i+1)).ge.abs(u(j,i+2)).and.
     #                abs(u(j,i+1)).gt.abs(u(j,i))) then
                     call vmove(t1,u(1,i),n)
                     call vmove(u(1,i),u(1,i+1),n)
                     call vmove(u(1,i+1),t1,n)
                  else if (abs(u(j,i+2)).ge.abs(u(j,i+1)).and.
     #                     abs(u(j,i+2)).gt.abs(u(j,i))) then
                     call vmove(t1,u(1,i),n)
                     call vmove(u(1,i),u(1,i+2),n)
                     call vmove(u(1,i+2),t1,n)
                  end if
c
c                 ----- annihilate this element in the second vector -----
c
                  a=-u(j,i+1)/u(j,i)
                  b=1.0d+00/sqrt(1.0d+00+a**2)
                  a=a*b
                  do 4 k=1,n
                     t2(k)=a*u(k,i)+b*u(k,i+1)
    4             continue
c
c                 ----- project the second vector out of the first -----
c
                  a=sdot(n,u(1,i),1,t2,1)
                  do 5 k=1,n
                     t1(k)=u(k,i)-a*t2(k)
    5             continue
c
c                 ----- normalize the first vector -----
c
                  a=sdot(n,t1,1,t1,1)
                  if (a.lt.fuzz) then
                     call lnkerr('error projecting in degrsp')
                  end if
                  a=1.0d+00/sqrt(a)
                  do 6 k=1,n
                     u(k,i)=a*t1(k)
                     u(k,i+1)=t2(k)
    6             continue
c
c                 ----- annihilate this element in the third vector -----
c
                  a=-u(j,i+2)/u(j,i)
                  b=1.0d+00/sqrt(1.0d+00+a**2)
                  a=a*b
                  do 7 k=1,n
                     t2(k)=a*u(k,i)+b*u(k,i+2)
    7             continue
c
c                 ----- project the third vector out of the first -----
c
                  a=sdot(n,u(1,i),1,t2,1)
                  do 8 k=1,n
                     t1(k)=u(k,i)-a*t2(k)
    8             continue
c
c                 ----- normalize the first vector -----
c
                  a=sdot(n,t1,1,t1,1)
                  if (a.lt.fuzz) then
                     call lnkerr('error projecting in degrsp')
                  end if
                  a=1.0d+00/sqrt(a)
                  do 9 k=1,n
                     u(k,i)=a*t1(k)
                     u(k,i+2)=t2(k)
    9             continue
c
c                 ----- increment the vector number by 1 so we fix the second
c                       pair of vectors
c
                  i=i+1
               end if
            end if
c
c           ----- degenerate pair or the second pair of a degenerate triplet
c
            do 10 j=1,n
               if (abs(u(j,i))+abs(u(j,i+1)).gt.small) go to 11
   10       continue
            call lnkerr('no value larger than small ????')
   11       continue
c
c           ----- switch the vectors so the one with the largest element is
c                 first
c
            if (abs(u(j,i)).lt.abs(u(j,i+1))) then
               call vmove(t1,u(1,i),n)
               call vmove(u(1,i),u(1,i+1),n)
               call vmove(u(1,i+1),t1,n)
            end if
c
c           ----- annihilate this coefficient in the second vector -----
c
            a=-u(j,i+1)/u(j,i)
            b=1.0d+00/sqrt(1.0d+00+a**2)
            a=a*b
            do 12 k=1,n
               t2(k)=a*u(k,i)+b*u(k,i+1)
   12       continue
c
c           ----- project the second vector out of the first -----
c
            a=sdot(n,u(1,i),1,t2,1)
            do 13 k=1,n
               t1(k)=u(k,i)-a*t2(k)
   13       continue
c
c           ----- and normalize the first vector -----
c
            a=sdot(n,t1,1,t1,1)
            if (a.lt.fuzz) call lnkerr('error projecting in degrsp')
            a=1.0d+00/sqrt(a)
            do 14 k=1,n
               u(k,i)=a*t1(k)
               u(k,i+1)=t2(k)
   14       continue
c
c           ----- increment the pointer to next vector -----
c
            i=i+1
         end if
         i=i+1
      if (i.lt.n) go to 1
c
c
      return
      end
