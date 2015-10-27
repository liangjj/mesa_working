*deck @(#)mcjac.f	1.1  11/30/90
      subroutine mcjac (dab,hess,xjk,nact,mrs,itran,thrsh)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcjac.f	1.1   11/30/90
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit real*8(a-h,o-z)
cc
cmp   extended dummy dab,hess,xjk
cc
      dimension dab(2),hess(2),xjk(2)
c
      if(itran.ne.0)go to 200
c
      nx=0
      do 120 i=1,nact
         den=-dab(i)
         if(abs(den).lt.thrsh)go to 115
         do 110 j=1,mrs
            hess(nx+j)=hess(nx+j)+den*xjk(j)
 110     continue
 115     continue
         nx=nx+mrs
 120  continue
c
      return
c
 200  continue
c
      do 220 j=1,mrs
         xint=xjk(j)
         if(abs(xint).lt.thrsh)go to 220
         nx=j
         do 210 i=1,nact
            hess(nx)=hess(nx)-xint*dab(i)
            nx=nx+mrs
 210     continue
 220  continue
c
      return
      end
