*deck @(#)sqtotr.f	5.1  11/6/94
      subroutine sqtotr(triang,square,num,nnp)
c***begin prologue     sqtotr
c***date written       850601  (yymmdd)
c***revision date      910729  (yymmdd)
c
c   29   july  1991    rlm at lanl
c      modifying so that if the sqare matrix is not symmetric
c      only the first bad element is printed and not the entire
c      matrix.
c***keywords           matrix, pack
c***author             saxe, paul (lanl)
c***source             @(#)sqtotr.f	5.1   11/6/94
c***purpose            packs a real symmetric matrix into a linear array.
c                      call sqtotr(triang,square,num,nnp)
c                        triang   output matrix, nnp words long.
c                        square   input matrix, dimensioned (num,num).
c                        nnp      num*(num+1)/2
c
c***references
c***routines called    abs
c***end prologue       sqtotr
      implicit integer (a-z)
c
      real*8 triang(nnp),square(num,num)
      logical notify
c
      common/io/inp,iout
c
      notify=.true.
      ij=0
      do 2 j=1,num
         do 1 i=1,j
            if (abs(square(i,j)-square(j,i)).gt.1.0d-05) then
               if (abs(square(i,j)-square(j,i))/
     #            max(abs(square(i,j)),abs(square(j,i)))
     #                 .gt.1.0d-05) then
               if(notify) then
                  write (iout,9) i,j,square(i,j),square(j,i)
               endif
               notify=.false.
    9          format (//,' ##### library: sqtotr, the square',
     #                 ' matrix is not symmetric:',2i4,2g18.9,//)
ctemp               call lnkerr('non-symmetric matrix')
               end if
            end if
c
            ij=ij+1
            triang(ij)=square(i,j)
    1    continue
    2 continue
c
      return
      end
