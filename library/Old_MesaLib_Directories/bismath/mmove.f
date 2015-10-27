*deck @(#)mmove.f
      subroutine mmove(a,b,n,m,ia,ib)
c***begin prologue     mmove
c***date written       960311  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix, move, transfer, copy
c***author             schneider, barry (nsf)
c***source             @(#)mmove.f
c***purpose            vectorized copy: a=b .
c***description
c                      call mmove(a,b,n,m,ia,ib)
c                        a        input matrix a(ia,*).
c                        b        output matrix b(ib,*).
c                        n,m      row and column indices.
c
c***references
c***routines called    (none)
c***end prologue       mmove
      real*8 a(ia,*), b(ib,*)
c
      do 1 i=1,n
         do 2 j=1,m
            b(i,j)=a(i,j)
 2       continue   
 1    continue
c
      return
      end
