*deck @(#)prthss.f	5.1  11/6/94
      subroutine prthss(buf,lbuf,tbuf,nmix)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)prthss.f	5.1   11/6/94
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
      dimension buf(lbuf),tbuf(*)
c
      common /io/ inp,iout
c
c
      write(iout,*)'  out-of-core hessian print  nmix ',nmix
c
      call iosys('rewind mcscf_hessian on rwf',0,0,0,' ')
      ntot=nmix*nmix
      lb=min(lbuf,ntot)
      lr=(lb/nmix)
      if(lr.lt.1) call lnkerr(' increase lbuf in prthss ')
      write(iout,*)'    read lb   ',lb
c
      write(iout,*)' finished read '
c
      jn=0
      in=lb
c
      is=0
      ie=0
c
      do 1 i=1,nmix,lr
         write(iout,*)'  row  ',i
         is=ie+1
         ie=is+lr
         ie=min(ie,nmix)
         len=ie-is+1
         lb=len*nmix
c
         call iosys('read real mcscf_hessian from rwf '//
     $        'without rewinding',lb,buf,0,' ')
c
         do 3 j=is,ie
            write(iout,*)'  row ',j
            write(iout,10)(buf(ix+jn),jn=1,nmix)
            ix=ix+nmix
 3       continue
c
 1    continue
c
 10   format(4(2x,f12.8))
      return
      end
