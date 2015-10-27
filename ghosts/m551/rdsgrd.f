*deck %W%  %G%
      subroutine rdsgrd(ahess,grad,lenb,nmix,nmx)
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
      dimension ahess(nmx,*),grad(nmix)
c
      data debug/.false./
      common /io/ inp,iout
c
c
      if(lenb.lt.nmix) then
         write(iout,*)'  lenb  nmix ',lenb,nmix
         call lnkerr('rdsgrd:increase buffer size(lenb)in mcaugh')
      endif
c
      call iosys('rewind mcscf_gradient on rwf',0,0,0,' ')
c
      call iosys('read real mcscf_gradient from rwf without rewinding',
     $     nmix,grad,0,' ')
c
      do 1 k=1,nmix
         ahess(k+1,1)=-grad(k)
 1    continue
c
      do 2 k=1,nmix
         ahess(1,k+1)=-grad(k)
 2    continue
c
      ahess(1,1)=0.0d+00
c
      if (debug) then
         write(iout,*) '  rdsgrd: ahess'
         call matout(ahess,nmx,nmx,nmx,nmx,iout)
      end if
c
      return
      end
