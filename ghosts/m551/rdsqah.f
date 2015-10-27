*deck %W%  %G%
      subroutine rdsqah(ahess,buf,lenb,nmix,nmx)
c
c***begin prologue
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
c***end prologue
c
      implicit real*8(a-h,o-z)
      logical debug
      dimension ahess(nmx,*),buf(nmix,*)
c
      data debug/.false./
      common /io/ inp,iout
c
c
      lr=(lenb/nmix)
      lr=min(lr,nmix)
c
      if(lr.lt.1) then
         write(iout,*)'  lenb  nmix ',lenb,nmix
         call lnkerr('rdsqah:increase buffer size(lenb)in mcaugh')
      endif
c
      call iosys('rewind mcscf_hessian on rwf',0,0,0,' ')
c
      ie=0
c
      do 1 i=1,nmix,lr
c
         lrr=min(lr,nmix-i+1)
         lb=lrr*nmix
         call iosys('read real mcscf_hessian from rwf '//
     $        'without rewinding',lb,buf,0,' ')
c
         is=ie+1
         ie=ie+lrr
         jx=0
         do 2 j=is,ie
            jx=jx+1
            do 3 k=1,nmix
               ahess(k+1,j+1)=buf(k,jx)
 3          continue
 2       continue
 1    continue
c
      if (debug) then
         write(iout,*) '  rdsqah: ahess'
         call matout(ahess,nmx,nmx,nmx,nmx,iout)
      end if
c
      return
      end
