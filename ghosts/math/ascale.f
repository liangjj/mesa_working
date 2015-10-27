*deck @(#)ascale.f	1.1  11/30/90
      subroutine ascale(n,s,a,b)
c***begin prologue     ascale
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, scale, scalar, multiply
c***author             martin, richard (lanl)
c***source             @(#)ascale.f	1.1   11/30/90
c***purpose            vectorized vector operation:
c                      b=s*a, where a,b are vectors, s a scalar.
c***description
c                      call ascale(n,s,a,b)
c                        n        vector length.
c                        s        scalar factor.
c                        a        input vector.
c                        b        output vector.
c
c***references
c***routines called    (none)
c***end prologue       ascale
      implicit real*8(a-h,o-z)
      dimension a(1),b(1)
c
c
      do 10 i=1,n
   10 b(i)=s*a(i)
c
      return
      end
