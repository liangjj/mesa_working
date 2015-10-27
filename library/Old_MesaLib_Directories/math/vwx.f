*deck @(#)vwx.f	5.1  11/6/94
      subroutine vwx(v,w,x,switch,n)
c***begin prologue     vwx
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, switch, add, subtract, negate, copy, move
c***author             saxe, paul (lanl)
c***source             @(#)vwx.f	5.1   11/6/94
c***purpose            vectorized operations: v=w+x, v=w-x, v=x, v=-x.
c***description
c                      call vwx(v,w,x,switch,n)
c                        v         output vector.
c                        w         input vector.
c                        x         input vector.
c                        switch    operation flag.
c                                  if  +1      0      -1        -2
c                                     v=w+x   v=x    v=w-x     v=-x
c                        n         vector length.
c
c***references
c***routines called    (none)
c***end prologue       vwx
cvax  extended dummy v,w,x
c
      real*8 v(n),w(n),x(n)
      integer switch
c
      if (switch.eq.0) then
         do 1 i=1,n
            v(i)=x(i)
    1    continue
      else if (switch.eq.1) then
         do 2 i=1,n
            v(i)=w(i)+x(i)
    2    continue
      else if (switch.eq.-1) then
         do 3 i=1,n
            v(i)=w(i)-x(i)
    3    continue
      else if (switch.eq.-2) then
         do 4 i=1,n
            v(i)=-x(i)
    4    continue
      else
         call lnkerr('unknown switch in vwx')
      end if
c
      return
      end
