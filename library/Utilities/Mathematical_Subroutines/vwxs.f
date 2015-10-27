*deck @(#)vwxs.f	5.1  11/6/94
      subroutine vwxs(v,w,x,scalar,switch,n)
c***begin prologue     vwxs
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, scalar, switch, add, subtract, negate, copy,
c                      move
c***author             saxe, paul (lanl)
c***source             @(#)vwxs.f	5.1   11/6/94
c***purpose            vectorized operations: v=w+x*s, v=w-x*s, v=x*s, v=-x*s,
c                      where v,w,x are vectors, s a scalar.
c***description
c                      call vwxs(v,w,x,s,switch,n)
c                        v         output vector.
c                        s         input scalar.
c                                  if    +1       0        -1          -2
c                                     v=w+x*s   v=x*s    v=w-x*s     v=-x*s
c                        n         vector length.
c
c***references
c***routines called    (none)
c***end prologue       vwxs
cvax  extended dummy v,w,x,y
c
      real*8 v(n),w(n),x(n),scalar
      integer switch
c
      if (switch.eq.0) then
         do 1 i=1,n
            v(i)=x(i)*scalar
    1    continue
      else if (switch.eq.1) then
         do 2 i=1,n
            v(i)=w(i)+x(i)*scalar
    2    continue
      else if (switch.eq.-1) then
         do 3 i=1,n
            v(i)=w(i)-x(i)*scalar
    3    continue
      else if (switch.eq.-2) then
         do 4 i=1,n
            v(i)=-x(i)*scalar
    4    continue
      else
         call lnkerr('unknown switch in vwxs')
      end if
c
      return
      end
