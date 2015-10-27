*deck @(#)bigues.f	1.1  11/30/90
      subroutine bigues(buf,diag,smlham,ipt,n,nsml,nlook,lenb,nmx)
c
c***begin prologue     bigues
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)bigues.f	1.1   11/30/90
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       bigues
c
      implicit real*8(a-h,o-z)
cmp   extended dummy buf,diag,smlham,ipt
      dimension buf(n,n),smlham(2),diag(2),ipt(2)
c
      common /io/ inp,iout
c
cc
      big=1.d+25
      small=1.d+20
cc
      do 20 i=1,nsml
         xmin=small
         imin=0
         do 10 j=1,nlook
            if(diag(j).gt.xmin) go to 10
            imin=j
            xmin=diag(j)
 10      continue
         if(imin.eq.0)go to 1000
         ipt(i)=imin
         smlham(i)=diag(imin)
         diag(imin)=big
 20   continue
cc
      do 30 i=1,nsml
         diag(ipt(i))=smlham(i)
 30   continue
cc
      do 40 i=1,nsml
         imin=ipt(i)
         do 50 j=i,nsml
            if(ipt(j).lt.imin) then
               imin=ipt(j)
               ipt(j)=ipt(i)
               ipt(i)=imin
            endif
 50      continue
 40   continue
cc
      imax=ipt(nsml)
cc
      lrr=lenb/nmx
      lr=min(lrr,imax)
      if(lrr.lt.1) call lnkerr(' buffer error in bigues of mcaugh ')
c
      lrd=lr*nmx
c
      call iosys('rewind mc_augh on mcscr',0,0,0,' ')
      call iosys('read real mc_augh from mcscr without rewinding',
     $     lrd,buf,0,' ')
c
      ls=0
      le=lr
c
      ix=0
c
      do 110 i=1,nsml
c
         if(ipt(i).gt.le) then
 100        continue
            left=imax-le
            lr=min(left,lrr)
            lrd=lr*nmx
            call iosys('read real mc_augh from mcscr without rewinding',
     $           lrd,buf,0,' ')
            ls=le
            le=le+lr
            if(ipt(i).gt.le) go to 100
         endif
c
         ii=ipt(i)-ls
         do 120 j=1,i
            ix=ix+1
            smlham(ix)=buf(ipt(j),ii)
 120     continue
 110  continue
cc
      return
cc
 1000 continue
      write(iout,1010)
 1010 format('  imin = 0 in guess   bug .. stop ')
      call lnkerr(' ')
      end
