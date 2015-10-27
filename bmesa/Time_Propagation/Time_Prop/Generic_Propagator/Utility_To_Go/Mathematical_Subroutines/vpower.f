*deck @(#)vpower.f	5.1  11/6/94
      subroutine vpower(v,w,ni,n)
c***begin prologue     vpower
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, power
c***author             saxe, paul (lanl)
c***source             @(#)vpower.f	5.1   11/6/94
c***purpose                                                       n
c                      vectorized vector to an integer power:  v=w  .
c***description
c                      call vpower(v,w,ni,n)
c                        v        output vector.
c                        w        input vector.
c                        ni       integer power to which elements of w are
c                                 raised.
c                        n        length of vector.
c
c***references
c***routines called    (none)
c***end prologue       vpower
      real*8 v(n),w(n)
      do 1 i=1,n
         v(i)=w(i)**ni
    1 continue
      return
      end
