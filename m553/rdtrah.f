*deck @(#)rdtrah.f	5.1  11/6/94
      subroutine rdtrah(ahess,buf,lenb,nmix,nmx)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)rdtrah.f	5.1   11/6/94
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
      implicit real*8 (a-h,o-z)
      logical debug
      dimension ahess(*),buf(nmix,*)
c
      parameter (debug=.false.)
      common /io/ inp,iout
c
c
      lr=(lenb/nmix)
      lr=min(lr,nmix)
c
      if(lr.lt.1) then
         write(iout,*)'rdtrah: lenb  nmix ',lenb,nmix
         call lnkerr('rdtrah:increase buffer size(lenb)in mcaugh')
      endif
c
      call iosys('rewind mcscf_hessian on rwf',0,0,0,' ')
c
      ie=0
      ix=1
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
            ix=ix+1
            jx=jx+1
            do 3 k=1,j
               ahess(ix+k)=buf(k,jx)
 3          continue
            ix=ix+j
 2       continue
 1    continue
c
      if (debug) then
         npp=nmx*(nmx+1)/2
         write(iout,9000) (ahess(i),i=1,npp)
 9000    format(/,'  rdtrah: ahess'/4(2x,f14.8))
      end if
c
      return
      end
