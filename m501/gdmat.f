*deck @(#)gdmat.f	5.1  11/6/94
      subroutine gdmat(d,c,nbf,nnp,mink,maxk)
c
c***module to form the density matrices
c
c paul saxe                  6 september 1984                 lanl
c
      implicit integer (a-z)
c
      real*8 d(nnp),c(nbf,nbf)
c
c
c
      ij=0
      call rzero(d,nnp)
      do 5 i=1,nbf
         do 4 j=1,i
            ij=ij+1
               do 1 k=mink,maxk
                  d(ij)=d(ij)+c(i,k)*c(j,k)
    1          continue
    4    continue
    5 continue
c
c
      return
      end
