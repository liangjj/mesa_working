*deck @(#)vexp.f	5.1  11/6/94
      subroutine vexp(v,w,n)
c***begin prologue     vexp
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, exponentiation
c***author             saxe, paul (lanl)
c***source             @(#)vexp.f	5.1   11/6/94
c***purpose            vectorized vector exponentiation:  v=exp(w) .
c***description
c                      call vexp(v,w,n)
c                        v        output vector of length (n).
c                        w        input vector of length (n).
c                        n        vector length.
c
c***references
c***routines called    (none)
c***end prologue       vexp
      real*8 v(n),w(n)
      do 1 i=1,n
         v(i)=exp(w(i))
    1 continue
      return
      end
