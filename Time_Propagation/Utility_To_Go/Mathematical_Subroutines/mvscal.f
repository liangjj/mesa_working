*deck @(#)mvscal.f	1.1  11/30/90
      subroutine mvscal(a,b,alpha,n)
c***begin prologue     mvscal
c***date written       921126  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           scale
c***author             schneider, barry (nsf)
c***source             @(#)mvscal.f
c***purpose            vectorized scale:  a=alpha*b
c***description
c                        a        output vector n.
c                        b        input matrix n*n.
c                        alpha    input scalar.
c
c***references
c***routines called    (none)
c***end prologue       mvscal
      real*8 a(n), b(n,n), alpha
      do 1 i=1,n
         a(i) = alpha*b(n,i)
    1 continue
      return
      end
