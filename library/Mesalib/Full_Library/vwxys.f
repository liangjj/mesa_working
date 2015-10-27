*deck @(#)vwxys.f	1.1  11/30/90
      subroutine vwxys(v,w,x,y,s,switch,n)
c***begin prologue     vwxys
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, switch, add, subtract, negate, multiply
c***author             saxe, paul (lanl)
c***source             @(#)vwxy.f	1.1   11/30/90
c***purpose            vectorized operations: v=w+s*x*y, v=w-s*x*y, v=s*x*y, 
c***                                          v=-s*x*y,
c                      where v,w,x, and y are vectors and s is a scalar.
c***description
c                      call vwxys(v,w,x,y,s,switch,n)
c                                  if    +1       0        
c                                     v=w+s*xy   v=s*xy    
c                        n         vector length.
c
c***references
c***routines called    (none)
c***end prologue       vwxys
c
      real*8 v(n), w(n), x(n), y(n), s
      integer switch
c
      if (switch.eq.0) then
         do 1 i=1,n
            v(i)=s*x(i)*y(i)
    1    continue
      else if (switch.eq.1) then
         do 2 i=1,n
            v(i)=w(i)+s*x(i)*y(i)
    2    continue
      else
         call lnkerr('unknown switch in vwxys')
      end if
c
      return
      end
