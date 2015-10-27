*deck %W%  %G%
      subroutine formah(ah,buf,lenb,grad,diag,nmix,nmx,smlld,shftd)
c
c***begin prologue     formah
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
c***end prologue       formah
c
      implicit real*8 (a-h,o-z)
      character*3 ians
      dimension grad(nmix),buf(nmix,*),ah(nmx,*),diag(nmx)
c
      common /io/ inp,iout
c
      call iosys('rewind mcscf_gradient on rwf',0,0,0,' ')
      call iosys('read real mcscf_gradient from rwf without rewinding',
     $     nmix,grad,0,' ')
c
      call iosys('does mc_augh exist on mcscr',0,0,0,ians)
c
      if(ians.eq.'no') then
         call iosys('create real mc_augh on mcscr',nmx*nmx,0,0,' ')
      endif
c
      call iosys('rewind mcscf_hessian on rwf',0,0,0,' ')
      call iosys('rewind mc_augh on mcscr',0,0,0,' ')
c
      lr=lenb/nmx
      lr=min(lr,nmx)
c
      if(lr.lt.1) then
         write(iout,*)' lenb must be > or = nmx    lenb nmx ',lenb,nmx
         call lnkerr(' increase buffer space in mcaugh')
      endif
c
      lwr=lr*nmx
      lrd=lr*nmix
c
      ah(1,1)=0.0d+00
      do 1 i=1,nmix
         ah(i+1,1)=-grad(i)
 1    continue
c
      if(lr.eq.1) then
         diag(1)=0.
         call iosys('write real mc_augh on mcscr without rewinding',
     $        nmx,ah,0,' ')
         lm=0
      else
         lm=lr-1
         lrr=lm*nmix
         call iosys('read real mcscf_hessian from rwf '//
     $        'without rewinding',lrr,buf,0,' ')
         ix=1
         do 2 i=1,lm
            ix=ix+1
            ah(1,ix)=-grad(i)
            do 3 j=1,nmix
               ah(j+1,ix)=buf(j,i)
 3          continue
 2       continue
         diag(1)=0.
         js=1
         lb=lm
         call bigdmp(ah(1,2),nmx,js,lb,smlld,shftd,diag(js+1))
         call iosys('write real mc_augh on mcscr without rewinding',
     $        lwr,ah,0,' ')
      endif
c
      left=nmix-lm
      if(left.gt.0) then
         jx=lm
         do 4 i=1,left,lr
            lb=min(lr,left-i+1)
            lrd=lb*nmix
            lwr=lb*nmx
            call iosys('read real mcscf_hessian from rwf '//
     $           'without rewinding',
     $           lrd,buf,0,' ')
c
            js=jx+1
c
            do 5 j=1,lb
               jx=jx+1
               ah(1,j)=-grad(jx)
               do 6 k=1,nmix
                  ah(k+1,j)=buf(k,j)
 6             continue
 5          continue
            call bigdmp(ah,nmx,js,lb,smlld,shftd,diag(js+1))
            call iosys('write real mc_augh to mcscr without rewinding',
     $           lwr,ah,0,' ')
 4       continue
      endif
c
c
      return
      end
