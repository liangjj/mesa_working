*deck %W%  %G%
      subroutine rdtgrd(ahess,grad,lenb,nmix,nmx)
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
      implicit real*8 (a-h,o-z)
      logical debug
      dimension ahess(*),grad(nmix)
c
      data debug/.false./
      common /io/ inp,iout
c
c
      if(lenb.lt.nmix) then
         write(iout,*)' rdtgrd: lenb  nmix ',lenb,nmix
         call lnkerr('rdtgrd:increase buffer size(lenb)in mcaugh')
      endif
c
      call iosys('rewind mcscf_gradient on rwf',0,0,0,' ')
c
      call iosys('read real mcscf_gradient from rwf without rewinding',
     $     nmix,grad,0,' ')
c
      ix=2
      do 1 k=1,nmix
         ahess(ix)=-grad(k)
         ix=ix+k+1
 1    continue
c
      ahess(1)=0.0d+00
c
      if (debug) then
         npp=nmx*(nmx+1)/2
         write(iout,9000) (ahess(i),i=1,npp)
 9000    format(/,'  rdtgrd: ahess '/4(2x,f14.8))
      end if
c
      return
      end
