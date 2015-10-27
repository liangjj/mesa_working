*deck @(#)sigmad.f	5.1  11/6/94
      subroutine sigmad(t,naxis,told)
c
      implicit integer (a-z)
c
      real*8 t(3,3,naxis),told(3,3,naxis),alpha,pi
      real*8 trans(3,3)
c
      pi=4*atan(1.0d+00)
      alpha=2*pi/(2*naxis)
c
      call rzero(trans,9)
      trans(1,1)=cos(alpha)
      trans(1,2)=sin(alpha)
      trans(2,1)=sin(alpha)
      trans(2,2)=-cos(alpha)
      trans(3,3)=1.0d+00
c
      do 1 i=1,naxis
         call ebc(t(1,1,i),trans,told(1,1,i),3,3,3)
    1 continue
c
      return
      end
