*deck %W%  %G%
      subroutine gethes(hess,lhess,buf,lenb,nf41)
c
c***begin prologue     gethes
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
c***end prologue       gethes
c
      implicit integer(a-z)
      real*8 hess(lhess),buf(lenb)
c
c
      call iosys('rewind mcscf_hessian on rwf',0,0,0,' ')
      call iosys('read real mcscf_hessian from rwf',lhess,hess,0,' ')
c
c
      return
      end
