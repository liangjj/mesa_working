      subroutine addrab(rab,buf,lbufso,tfile,noc,nob,nao)
C
C***Begin prologue     addrab
C***Date written       871022   (yymmdd)
C***Revision date      yymmdd   (yymmdd)
C
C***Keywords
C***Author             Lengsfield, Byron (BRL)
C***Source             %W%   %G%
C
C***Purpose
C
C***Description
C
C***References
C
C***Routines called    (none)
C
C***End prologue       addrab
C
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
      call iosys('rewind '//file//' on rwf',0,0,0,0)
c
      nni=naa
      ix=1
c
      do 10 i=1,npass
c
         nleft=min(lpass,nni)
         nni=nni-nleft
c
         nr=nleft*noi
         call iosys('read real '//file//' from rwf without rewinding',
     $        lbufso,buf,0,0)
c
         call vadd(rab(ix),rab(ix),buf,nr)
         ix=ix+nr
c
 10   continue
c
      return
      end
