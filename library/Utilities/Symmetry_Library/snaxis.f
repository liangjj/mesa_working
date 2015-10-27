*deck @(#)snaxis.f	5.1  11/6/94
      subroutine snaxis(t,naxis)
c
c     the functions transform as a clockwise rotation;(2*pi*m/naxis).
c     followed by a reflection perpendicualr to the plane.
c
      implicit integer (a-z)
c
      real*8 t(3,3,0:*),alpha,pi
c
      pi=atan(1.0d+00)*4
      alpha=2.0d+00*pi/naxis
c
      call rzero(t,3*3*naxis)
      do 1 i=0,naxis-1
         t(1,1,i)=-cos(i*alpha)
         t(1,2,i)=-sin(i*alpha)
         t(2,1,i)=+sin(i*alpha)
         t(2,2,i)=-cos(i*alpha)
         t(3,3,i)= 1.0d+00
    1 continue
c
c
      return
      end
