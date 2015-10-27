*deck @(#)fixdm.f	5.1  11/6/94
      subroutine fixdm(h,g,norbs,nnp)
c
      implicit integer (a-z)
c
      real*8 h(nnp),g(nnp,nnp)
c
c     ----- eliminate all factors of two and four from replication
c
      do 2 i=1,norbs
         do 1 j=1,i-1
            ij=i*(i-1)/2+j
            h(ij)=h(ij)*0.5d+00
    1    continue
    2 continue
c
      do 4 j=1,nnp
         do 3 i=1,nnp
            g(i,j)=g(i,j)*0.25d+00
    3    continue
    4 continue
c
      ii=0
      do 7 i=1,norbs
         ii=ii+i
         do 5 j=1,nnp
            g(j,ii)=g(j,ii)*2.0d+00
    5    continue
         do 6 j=1,nnp
            g(ii,j)=g(ii,j)*2.0d+00
    6    continue
    7 continue
c
c     ----- scale g by 1/2  -----
c
      do 12 j=1,nnp
         do 11 i=1,nnp
            g(i,j)=g(i,j)*0.5d+00
   11    continue
   12 continue
c
c     ----- g(ik,kj) = h(ik,kj) - h(ij) -----
c
      ij=0
      do 10 i=1,norbs
         do 9 j=1,i
            ij=ij+1
            do 8 k=1,norbs
               ii=max(i,k)
               ik=ii*(ii-1)/2+min(i,k)
               kk=max(j,k)
               kj=kk*(kk-1)/2+min(j,k)
               if (i.eq.j) then
                  if (i.eq.k.and.k.eq.j) then
                     g(ik,kj)=g(ik,kj)-h(ij)*0.50d+00
                  else if (i.eq.k.or.k.eq.j) then
                     g(ik,kj)=g(ik,kj)-h(ij)*0.25d+00
                  else
                     g(ik,kj)=g(ik,kj)-h(ij)*0.125d+00
                  end if
               else
                  if (i.eq.k.or.k.eq.j) then
                     g(ik,kj)=g(ik,kj)-h(ij)*0.25d+00
                     g(kj,ik)=g(kj,ik)-h(ij)*0.25d+00
                  else
                     g(ik,kj)=g(ik,kj)-h(ij)*0.125d+00
                     g(kj,ik)=g(kj,ik)-h(ij)*0.125d+00
                  end if
               end if
cps               if (i.ne.j) then
cps                  g(kj,ik)=g(kj,ik)-h(ij)
cps               end if
    8       continue
    9    continue
   10 continue
c
c
      return
      end
