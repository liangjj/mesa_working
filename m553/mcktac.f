*deck @(#)mcktac.f	5.1  11/6/94
      subroutine mcktac(dab,hess,xjk,nact,nr,ns,itran,thrsh,fact)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcktac.f	5.1   11/6/94
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
      mrs=nr*ns
      if(itran.ne.0)go to 200
c
      nx=0
      do 120 i=1,nact
         den=dab(i)*fact
         if(abs(den).lt.thrsh)go to 115
         nxe=0
         do 111 k=1,ns
            kx=k
            nxs=nxe+1
            nxe=nxe+nr
            do 110 j=nxs,nxe
               hess(nx+j)=hess(nx+j)+den*xjk(kx)
               kx=kx+nr
 110        continue
 111     continue
 115     continue
         nx=nx+mrs
 120  continue
c
      return
c
 200  continue
c
      mx=0
      do 230 k=1,ns
         kx=k
         do 220 j=1,nr
            xint=xjk(kx)*fact
            kx=kx+nr
            mx=mx+1
            if(abs(xint).lt.thrsh)go to 220
            nx=mx
            do 210 i=1,nact
               hess(nx)=hess(nx)+xint*dab(i)
               nx=nx+mrs
 210        continue
 220     continue
 230  continue
c
      return
      end
