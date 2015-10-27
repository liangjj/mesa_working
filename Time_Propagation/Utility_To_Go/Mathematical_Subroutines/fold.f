*deck @(#)fold.f	5.1  11/6/94
      subroutine fold(triang,square,num,nnp)
c***begin prologue     fold
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix, pack
c***author             saxe, paul (lanl)
c***source             @(#)fold.f	5.1   11/6/94
c***purpose            packs a real symmetric matrix into a linear array.
c                      call fold(triang,square,num,nnp)
c                        triang   output matrix, nnp words long.
c                        square   input matrix, dimensioned (num,num).
c                        nnp      num*(num+1)/2
c
c***references
c***routines called    abs
c***end prologue       fold
      implicit integer (a-z)
c
      real*8 triang(nnp),square(num,num)
c
      common/io/inp,iout
c
      ij=0
      do 2 j=1,num
         do 1 i=1,j-1
            ij=ij+1
            triang(ij)=(square(i,j)+square(j,i))/2.0d+00
    1    continue
         ij=ij+1
         triang(ij)=square(j,j)
    2 continue
c
      return
      end
