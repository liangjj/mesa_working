*deck rdpoly.f
c***begin prologue     rdpoly
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            lobatto information.
c***                   
c***references         
c
c***routines called    
c***end prologue       rdpoly
      subroutine rdpoly(px,pxply,xloc,pxloc,ngrid,npt,ngot)
      implicit integer (a-z)
      real*8 x, xply
      character*2 itoc, ig
      dimension npt(ngrid), ngot(2), xloc(ngrid), pxloc(ngrid,ngrid)
      common/io/inp, iout
      pointer (px,x(1))
      pointer (pxply,xply(1))
c
c     compute memory pointers and needed memory
c
      cnti=1
      cntji=1
      do 10 grdi=1,ngrid
         xloc(grdi)=cnti
         xi=cnti
         xwti=q+npt(grdi)
         cnti=xwti+npt(grdi)
         do 20 grdj=1,ngrid
            pxloc(grdj,grdi)=cntji
            pji=cntji
            dpji=pji+npt(grdi)*npt(grdj)
            ddpji=dpji+npt(grdi)*npt(grdj)
            cntji=ddpji+npt(grdi)*npt(grdj)
 20      continue
 10   continue
c
c     get the memory
c 
      cnti=wpadti(cnti)
      cntji=wpadti(cntji)
      call memory(cnti,px,ngot(1),'x',0)
      call memory(cntji,pxply,ngot(2),'pji',0)
c
c     read in the data
c
      do 30 grdi=1,ngrid
         ig=itoc(grdi)
         xi=xloc(grdi)
         xwti=xi+npt(grdi)
         call iosys('read real "points grid = '//ig//'" from bec',
     1               npt(grdi),x(xi),0,' ')
         call iosys('read real "weights grid = '//ig//'" from bec',
     1               npt(grd),x(xwti),0,' ')
         do 40 grdj=1,ngrid
            jg=itoc(grdj) 
            pji=pxloc(grdj,grdi)
            dpji=pji+npt(grdi)*npt(grdj)
            ddpji=dpji+npt(grdi)*npt(grdj)
            call iosys('read real "pji set = '//ig//' grid = '//jg
     1                  //'" from bec',npt(grdi)*npt(grdj),xply(pji),
     2                  0,' ')                  
            call iosys('read real "dpji set = '//ig//' grid = '//jg
     1                  //'" from bec',npt(grdi)*npt(grdj),xply(dpji),
     2                  0,' ')                  
            call iosys('read real "ddpji set = '//ig//' grid = '//jg
     1                  //'" from bec',npt(grdi)*npt(grdj),xply(ddpji),
     2                  0,' ')
 40      continue   
 30   continue   
      return
      end       


