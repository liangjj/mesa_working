*deck %W%  %G%
      subroutine bighvc(buf,c,b,nmx,nvec,lenb)
c
c***begin prologue     bighvc
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       bighvc
c
c
      implicit real*8 (a-h,o-z)
      dimension b(nmx,nvec),c(nmx,nvec),buf(nmx,*)
c
      lr=lenb/nmx
      lr=min(lr,nmx)
c
      if(lr.lt.1) call lnkerr(' increase lenb for bighvc in mcaugh')
c
      call iosys('rewind mc_augh on mcscr',0,0,0,' ')
c
      ix=1
      do 1 i=1,nmx,lr
c
         lb=min(lr,nmx-i+1)
         lrd=lb*nmx
c
         call iosys(' read real mc_augh from mcscr without rewinding',
     $        lrd,buf,0,' ')
ccc
c      b(k,l)=buf(j,k)*c(j,l)
ccc
         call mxma(buf,nmx,1,c,1,nmx,b(ix,1),1,nmx,lb,nmx,nvec)
c
         ix=ix+lb
c
 1    continue
c
      return
      end
