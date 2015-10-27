*deck @(#)trtosq.f	5.1  11/6/94
      subroutine trtosq(square,triang,num,nnp)
c***begin prologue     trtosq
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix, unpack
c***author             saxe, paul (lanl)
c***source             @(#)trtosq.f	5.1   11/6/94
c***purpose            unpacks a triangularly packed matrix into a square.
c***description
c                      call trtosq(square,triang,num,nnp)
c                        square   output matrix, dimensioned (num,num).
c                        triang   input matrix, nnp words long.
c                        nnp      num*(num+1)/2
c
c***references
c***routines called    (none)
c***end prologue       trtosq
      implicit integer (a-z)
c
      real*8 triang(nnp),square(num,num)
c
      ij=0
      do 2 j=1,num
         do 1 i=1,j
            ij=ij+1
            square(i,j)=triang(ij)
            square(j,i)=triang(ij)
    1    continue
    2 continue
c
c
      return
      end
