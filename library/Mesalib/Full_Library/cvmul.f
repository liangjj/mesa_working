*deck @(#)cvmul.f	1.1  11/30/90
      subroutine cvmul(v,w,xr,xc,n,xtype)
c***begin prologue     cvmul
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, multiply
c***author             saxe, paul (lanl)
c***source             @(#)vmul.f	1.1   11/30/90
c***purpose            vectorized vector multiply:  v=w*x .
c***description
c                      call cvmul(v,w,xr,xc,n)
c                        v        complex output vector of length n.
c                        w        complex input vector of length n.
c                        xr       real input vector of length n.
c                        xc       complex input vector of length n
c                        n        length of vector.
c
c***references
c***routines called    (none)
c***end prologue       vmul
      complex*16 v(n),w(n),xc(n)
      real*8 xr(n)
      character*(*) xtype
      if (xtype.eq.'real') then
          do 1 i=1,n
             v(i)=w(i)*xr(i)
    1     continue
      elseif(xtype.eq.'complex') then
	  do 2 i=1,n
	     v(i)=w(i)*xc(i)
    2     continue
      endif 
      return
      end
