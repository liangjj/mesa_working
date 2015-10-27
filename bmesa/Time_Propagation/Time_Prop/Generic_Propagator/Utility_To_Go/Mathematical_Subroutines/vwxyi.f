*deck @(#)vwxyi.f	5.1  11/6/94
      subroutine vwxyi(v,w,x,y,index,switch,n)
c***begin prologue     vwxyi
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, switch, add, subtract, negate, multiply, gather
c***author             saxe, paul (lanl)
c***source             @(#)vwxyi.f	5.1   11/6/94
c***purpose            vectorized operations: v=w+x*y, v=w-x*y, v=x*y, v=-x*y,
c                      where v,w,x, and y are vectors.  the elements of y are
c                      gathered via a pointer array.
c***description
c                        w         input vector.
c                        x         input vector.
c                        y         input vector.
c                        index     array pointing to the elements of y.
c                                  if    +1       0        -1          -2
c                                     v=w+xy   v=xy    v=w-xy     v=-xy
c***routines called    (none)
c***end prologue       vwxyi
cvax  extended dummy v,w,x,y,index
c
      real*8 v(n),w(n),x(n),y(*)
      integer switch,index(n)
c
      if (switch.eq.0) then
         do 1 i=1,n
            v(i)=x(i)*y(index(i))
    1    continue
      else if (switch.eq.1) then
         do 2 i=1,n
            v(i)=w(i)+x(i)*y(index(i))
    2    continue
      else if (switch.eq.-1) then
         do 3 i=1,n
            v(i)=w(i)-x(i)*y(index(i))
    3    continue
      else if (switch.eq.-2) then
         do 4 i=1,n
            v(i)=-x(i)*y(index(i))
    4    continue
      else
         call lnkerr('unknown switch in vwxyi')
      end if
c
      return
      end
