*deck @(#)csmul.f	1.1  11/30/90
      subroutine csmul(v,w,s,n)
c***begin prologue     csmul
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, scalar, multiply.
c***author             saxe, paul (lanl)
c***source             @(#)csmul.f	1.1   11/30/90
c***purpose            vectorized scalar multiply v=w*s, where v,w are vectors
c                      and s is scalar. variables are all complex.
c***description
c                      call csmul(v,w,s,n)
c                        v        output vector dimensioned (n).
c                        w        input vector dimensioned (n).
c                        s        scalar.
c
c***references
c***routines called    (none)
c***end prologue       csmul
      complex*16 v(n),w(n),s
      do 1 i=1,n
         v(i)=w(i)*s
    1 continue
      return
      end
