*deck @(#)rdiagh.f	1.1  11/30/90
      subroutine rdiagh(diag,buf,lenb,nmix)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)rdiagh.f	1.1   11/30/90
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
      dimension diag(*),buf(nmix,*)
c
      common /io/ inp,iout
c
c
      lr=(lenb/nmix)
      lr=min(lr,nmix)
c
      if(lr.lt.1) then
         write(iout,*)'  lenb  nmix ',lenb,nmix
         call lnkerr('rdiagh:increase buffer size(lenb)in mcaugh')
      endif
c
      call iosys('rewind mcscf_hessian on rwf',0,0,0,' ')
c
      ie=0
      ix=0
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
            diag(j)=buf(j,jx)
 2       continue
c
 1    continue
c
      return
      end
