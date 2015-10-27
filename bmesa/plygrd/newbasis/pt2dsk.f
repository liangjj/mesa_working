*deck pt2dsk.f
c***begin prologue     pt2dsk
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            disk io for lobatto points and weights.
c***                   
c***references         
c
c***routines called    
c***end prologue       pt2dsk
      subroutine pt2dsk(ppt,ngrid,npt,phrse,prn)
      implicit integer (a-z)
      real*8 pt
      dimension npt(ngrid)
      character*(*) phrse
      character*2 itoc, ig
      character*80 title
      common/io/inp, iout
      pointer (ppt,pt(1))
      title=phrse//'"'
      write(iout,1) title
      begin=1
      do 10 grd=1,ngrid
         q=begin
         wt=q+npt(grd)
         begin=wt+npt(grd)
         ig=itoc(grd)
         title=phrse//' grid = '//ig//' no. points"'
         call iosys('write real '//title//' to bec',1,npt(grd),0,' ')
         title=phrse//' grid = '//ig//' points"'
         call iosys('write real '//title//' to bec',npt(grd),
     1               pt(q),0,' ')
         if(prn) then
            call prntrm(title,pt(q),npt(grd),1,npt(grd),1,iout)
         endif
         title=phrse//' grid = '//ig//' weights"'
         call iosys('write real '//title//' to bec',npt(grd),
     1               pt(wt),0,' ')
         if(prn) then
            call prntrm(title,pt(wt),npt(grd),1,npt(grd),1,iout)
         endif
 10   continue   
      return
 1    format(/,10x,a32)
 2    format(/,10x,'grid              = ',i2,/10x,
     1             'number of points  = ',i3)
 2    format(/,10x,'grid                      = ',i2,/10x,
     1             'modified number of points = ',i3)
      end       
