*deck @(#)fmints.f	5.1  11/6/94
      subroutine fmints(h,g,norbs,nnp)
c
      implicit integer (a-z)
c
      real*8 h(nnp),g(nnp,nnp)
c
c     ----- scale g by 1/2  -----
c
      do 2 j=1,nnp
         do 1 i=1,nnp
            g(i,j)=g(i,j)*0.5d+00
    1    continue
    2 continue
c
c     ----- form h = h(ij)-1/2 sum(k) [ik;kj]   -----
c
      ij=0
      do 7 i=1,norbs
         do 6 j=1,i
            ij=ij+1
            do 5 k=1,norbs
               ii=max(i,k)
               ik=ii*(ii-1)/2+min(i,k)
               kk=max(j,k)
               kj=kk*(kk-1)/2+min(j,k)
               h(ij)=h(ij)-g(ik,kj)
    5       continue
    6    continue
    7 continue
c
c
      return
      end
