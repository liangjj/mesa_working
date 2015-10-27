*deck ply2dsk.f
c***begin prologue     ply2dsk
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            disk io for lobatto polynomials.
c***                   
c***references         
c
c***routines called    
c***end prologue       ply2dsk
      subroutine ply2dsk(ppoly,ngrid,npt,phrse,prn)
      implicit integer (a-z)
      real*8 poly
      dimension npt(ngrid)
      dimension title(2)
      character*(*) phrse
      character*2 itoc, ig, jg
      character*80 title
      common/io/inp, iout
      pointer (ppoly,poly(1))
      cntply=1
      do 10 grdi=1,ngrid
         ig=itoc(grdi)
         title(1)=phrse//' set = '//itoc(ig)  
         do 20 grdj=1,ngrid
            jg=itoc(grdj)
            p=cntply
            dp=p+npt(grdi)*npt(grdj)
            ddp=dp+npt(grdi)*npt(grdj)
            cntply=ddp+npt(grdi)*npt(grdj)
            title(2)=title(1)//' grid = '//jg//' pji"'
            call iosys('write real '//title(2)//' to bec',
     1                  npt(grdi)*npt(grdj),poly(p),0,' ')
            if(prn) then
               call prntrm(title(2),poly(p),npt(grdj),npt(grdi),
     1                     npt(grdj),npt(grdi),iout)
            endif
            title(2)=title(1)//' grid = '//jg//' dpji"'
            call iosys('write real '//title(2)//' to bec',
     1                  npt(grdi)*npt(grdj),poly(dp),0,' ')
            if(prn) then
               call prntrm(title(2),poly(dp),npt(grdj),npt(grdi),
     1                     npt(grdj),npt(grdi),iout)
            endif
            title(2)=title(1)//' grid = '//jg//' ddpji"'
            call iosys('write real '//title(2)//' to bec',
     1                  npt(grdi)*npt(grdj),poly(ddp),0,' ')
            if(prn) then
               call prntrm(title(2),poly(ddp),npt(grdj),npt(grdi),
     1                     npt(grdj),npt(grdi),iout)
            endif
 20      continue   
 10   continue   
      return
      end       


