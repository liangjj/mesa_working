*deck @(#)smul.f	5.1  11/6/94
      subroutine smul(v,w,s,n)
c***begin prologue     smul
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, scalar, multiply.
c***author             saxe, paul (lanl)
c***source             @(#)smul.f	5.1   11/6/94
c***purpose            vectorized scalar multiply v=w*s, where v,w are vectors
c                      and s is scalar.
c***description
c                      call smul(v,w,s,n)
c                        v        output vector dimensioned (n).
c                        w        input vector dimensioned (n).
c                        s        scalar.
c
c***references
c***routines called    (none)
c***end prologue       smul
      real*8 v(n),w(n),s
      do 1 i=1,n
         v(i)=w(i)*s
    1 continue
      return
      end
