*deck @(#)sadd.f	5.1  11/6/94
      subroutine sadd(v,w,s,n)
c***begin prologue     sadd
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, scalar, add
c***author             saxe, paul (lanl)
c***source             @(#)sadd.f	5.1   11/6/94
c***purpose            vectorized scalar addition to vector:  v=w+s,
c                      where v,w are vectors and s a scalar.
c***description
c                      call sadd(v,w,s,n)
c                        v        output vector.
c                        w        input vector.
c                        s        scalar.
c
c***references
c***routines called    (none)
c***end prologue       sadd
      real*8 v(n),w(n),s
      do 1 i=1,n
         v(i)=w(i)+s
    1 continue
      return
      end
