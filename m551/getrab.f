*deck @(#)getrab.f	1.1  11/30/90
      subroutine getrab(rab,buf,lbufso,tfile,noc,nob,nao)
c
c***begin prologue     getrab
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)getrab.f	1.1   11/30/90
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       getrab
c
      implicit integer(a-z)
      real*8 buf(lbufso),rab(*)
c
      character*(*) tfile,file*16
c
      file=tfile
c
      naa=nao*(nao+1)/2
      noi=noc*nob
      lpass=lbufso/noi
      npass=(naa-1)/lpass+1
c
      call iosys('rewind '//file//' on rwf',0,0,0,' ')
c
      nni=naa
      ix=1
      do 10 i=1,npass
         nleft=min(lpass,nni)
         nni=nni-nleft
         nr=nleft*noi
         call iosys('read real '//file//' from rwf without rewinding',
     $        lbufso,buf,0,' ')
         call vmove(rab(ix),buf,nr)
         ix=ix+nr
 10   continue
c
      return
      end
