*deck @(#)mscale.f	1.1  11/30/90
      subroutine mscale(alpha,a,n,m)
c***begin prologue     mscale
c***date written       921126  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix, scale
c***author             schneider, barry (nsf)
c***source             @(#)mscale.f	1.1   11/26/92
c***purpose            vectorized matrix scale:  a=alpha*a
c***description
c                        a        input matrix n*m.
c                        alpha    input scalar.
c
c***references
c***routines called    (none)
c***end prologue       mscale
      real*8 a(n,m), alpha
      do 1 i=1,n
         do 2 j=1,m
            a(i,j) = alpha*a(i,j)
    2    continue
    1 continue
      return
      end
