*deck @(#)gdmat.f	5.1  11/6/94
      subroutine gdmat (d,c,num,newmo,nnp,mink,maxk)
c
c***module to form the density matrices
c
c paul saxe                  6 september 1984                 lanl
c
      implicit integer(a-z)
c
      real*8 d(nnp), c(num,newmo)
c
      ij=0
      call rzero (d,nnp)
      do 30 i=1,num
      do 20 j=1,i
      ij=ij+1
      do 10 k=mink,maxk
      d(ij)=d(ij)+c(i,k)*c(j,k)
   10 continue
   20 continue
   30 continue
c
      return
      end
