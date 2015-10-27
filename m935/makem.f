*deck @(#)makem.f	5.1  11/6/94
       subroutine makem(thess,r,liter,mvec,iter)
       implicit real*8(a-h,o-z)
       dimension thess(*),r(iter,mvec)
c
       if(liter.ne.0) then
          ix=1
          do 1 i=1,mvec
             call scopy(liter,r(1,i),1,thess(ix),1)
             ix=ix+liter
             call scopy(i,r(liter+1,i),1,thess(ix),1)
             ix=ix+i
  1       continue
       else
          ix=1
          do 2 i=1,mvec
             call scopy(i,r(1,i),1,thess(ix),1)
             ix=ix+i
  2       continue
       end if
c
       return
       end
