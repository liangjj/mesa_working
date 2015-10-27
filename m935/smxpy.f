*deck @(#)smxpy.f	5.1  11/6/94
      subroutine smxpy(n,y,nvec,ld,s,xm)
      implicit real*8(a-h,o-z)
      real*8 y(*),s(*),xm(ld,*)
c
      do 1 i=1,nvec
         call saxpy(n,s(i),xm(1,i),1,y,1)
  1   continue
c
      return
      end
