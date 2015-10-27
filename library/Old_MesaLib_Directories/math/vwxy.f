*deck @(#)vwxy.f	5.1  11/6/94
      subroutine vwxy(v,w,x,y,switch,n)
c***begin prologue     vwxy
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, switch, add, subtract, negate, multiply
c***author             saxe, paul (lanl)
c***source             @(#)vwxy.f	5.1   11/6/94
c***purpose            vectorized operations: v=w+x*y, v=w-x*y, v=x*y, v=-x*y,
c                      where v,w,x, and y are vectors.
c***description
c                      call vwxy(v,w,x,y,switch,n)
c                                  if    +1       0        -1          -2
c                                     v=w+xy   v=xy    v=w-xy     v=-xy
c                        n         vector length.
c
c***references
c***routines called    (none)
c***end prologue       vwxy
cvax  extended dummy v,w,x,y
c
      real*8 v(n),w(n),x(n),y(n)
      integer switch
c
      if (switch.eq.0) then
         do 1 i=1,n
            v(i)=x(i)*y(i)
    1    continue
      else if (switch.eq.1) then
         do 2 i=1,n
            v(i)=w(i)+x(i)*y(i)
    2    continue
      else if (switch.eq.-1) then
         do 3 i=1,n
            v(i)=w(i)-x(i)*y(i)
    3    continue
      else if (switch.eq.-2) then
         do 4 i=1,n
            v(i)=-x(i)*y(i)
    4    continue
      else
         call lnkerr('unknown switch in vwxy')
      end if
c
      return
      end
