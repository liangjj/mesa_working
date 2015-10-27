*deck @(#)enforc.f	3.1  11/20/92
      subroutine enforc(n,nnp,s,smhalf,lambda,u,t1,t2)
c
c***begin prologue     enforc
c***date written       860813  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           degenerate eigenvalues
c***author             saxe, paul (lanl)
c***source             @(#)enforc.f	3.1   11/20/92
c***purpose            enforcing symmetry in the scf vector by
c                      rotation of degenerate eigenvectors to sensible
c                      arrangement.
c***description
c                      call enforc(n,nnp,s,lambda,u,t1,t2)
c
c      n         integer
c                the matrix dimension
c
c      nnp       integer
c                the size of the packed matrix [n*(n+1)/2]
c
c      s         real (nnp)
c                the overlap matrix, packed.
c
c      smhalf    real (nnp)
c                s**-1/2 in packed form.
c
c      lambda    real (n)
c                the eigenvalues.
c
c      u         real (n,n)
c                the array for eigenvectors.
c
c      t1, t2    real (n,n)
c                scratch arrays.
c
c      this routine currently handles doubly and triply degenerate eigenvalues
c
c***references
c***routines called
c***end prologue       enforc
c
      implicit integer (a-z)
c
      real*8 s(nnp),u(n,n),lambda(n),t1(n),t2(n),smhalf(nnp)
      real*8 fuzz,reldif,a,b,small,sdot
c
      parameter (fuzz=1.0d-06,small=1.0d-01)
c
c     ----- statement functions -----
c
      reldif(a,b)=abs((a-b)/(a+b))
c
c     ----- transform to an orthonormal basis, namely s**-(1/2) . u
c
      call trtosq(t1,s,n,nnp)
      call ebc(t2,t1,u,n,n,n)
      call trtosq(t1,smhalf,n,nnp)
      call ebc(u,t1,t2,n,n,n)
c
c----------------------------------------------------------------------
c   rotate eigenvectors so that they are unique. strategy is to
c   subtract the first of a degenerate pair from the second with a
c   coefficient such that the first largish element of the second
c   eigenvector is annihilated. the second vector is then projected from
c   the first to hopefully leave two nice vectors.
c-----------------------------------------------------------------------
c
      i=1
    1 continue
         if (reldif(lambda(i),lambda(i+1)).lt.fuzz) then
            if (i.lt.n-1) then
               if (reldif(lambda(i+1),lambda(i+2)).lt.fuzz) then
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
     $                abs(u(j,i+1)).gt.abs(u(j,i))) then
                     call vmove(t1,u(1,i),n)
                     call vmove(u(1,i),u(1,i+1),n)
                     call vmove(u(1,i+1),t1,n)
                  else if (abs(u(j,i+2)).ge.abs(u(j,i+1)).and.
     $                     abs(u(j,i+2)).gt.abs(u(j,i))) then
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
                     u(k,i+1)=a*u(k,i)+b*u(k,i+1)
    4             continue
c
c                 ----- project the second vector out of the first -----
c
                  a=sdot(n,u(1,i),1,u(1,i+1),1)
                  do 5 k=1,n
                     u(k,i)=u(k,i)-a*u(k,i+1)
    5             continue
c
c                 ----- normalize the first vector -----
c
                  a=sdot(n,u(1,i),1,u(1,i),1)
                  if (a.lt.fuzz) call lnkerr('normalization error')
                  a=1.0d+00/sqrt(a)
                  do 6 k=1,n
                     u(k,i)=u(k,i)*a
    6             continue
c
c                 ----- annihilate this element in the third vector -----
c
                  a=-u(j,i+2)/u(j,i)
                  b=1.0d+00/sqrt(1.0d+00+a**2)
                  a=a*b
                  do 7 k=1,n
                     u(k,i+2)=a*u(k,i)+b*u(k,i+2)
    7             continue
c
c                 ----- project the third vector out of the first -----
c
                  a=sdot(n,u(1,i),1,u(1,i+2),1)
                  do 8 k=1,n
                     u(k,i)=u(k,i)-a*u(k,i+2)
    8             continue
c
c                 ----- normalize the first vector -----
c
                  a=sdot(n,u(1,i),1,u(1,i),1)
                  if (a.lt.fuzz) call lnkerr('normalization error')
                  a=1.0d+00/sqrt(a)
                  do 9 k=1,n
                     u(k,i)=u(k,i)*a
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
               u(k,i+1)=a*u(k,i)+b*u(k,i+1)
   12       continue
c
c           ----- project the second vector out of the first -----
c
            a=sdot(n,u(1,i),1,u(1,i+1),1)
            do 13 k=1,n
               u(k,i)=u(k,i)-a*u(k,i+1)
   13       continue
c
c           ----- normalize the first vector -----
c
            a=sdot(n,u(1,i),1,u(1,i),1)
            if (a.lt.fuzz) call lnkerr('normalization error')
            a=1.0d+00/sqrt(a)
            do 14 k=1,n
               u(k,i)=u(k,i)*a
   14       continue
c
c           ----- increment the pointer to next vector -----
c
            i=i+1
         end if
         i=i+1
      if (i.lt.n) go to 1
c
c     ----- transform back to the non-orthogonal basis -----
c
      call trtosq(t1,smhalf,n,nnp)
      call ebc(t2,t1,u,n,n,n)
      call vmove(u,t2,n**2)
c
c     ----- check orthonormality of result -----
c
      call trtosq(t1,s,n,nnp)
      call ebc(t2,t1,u,n,n,n)
      call ebtc(t1,u,t2,n,n,n)
c
      ij=0
      do 32 i=1,n
         do 31 j=1,n
            ij=ij+1
            if (i.ne.j) then
               if (abs(t1(ij)).gt.fuzz) call lnkerr('orthogonality')
            else
               if (abs(1.0d+00-t1(ij)).gt.fuzz) call lnkerr('norm')
            end if
   31    continue
   32 continue
c
c
      return
      end
