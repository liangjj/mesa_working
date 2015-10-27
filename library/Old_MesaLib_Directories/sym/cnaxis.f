*deck @(#)cnaxis.f	5.1  11/6/94
      subroutine cnaxis(t,naxis,axis)
c
c     the coordinates transform as a counter-clockwise rotation;(2*pi*m/naxis).
c     axis denotes whether the rotation is to be about x, y, or z(default).
c
      implicit integer (a-z)
c
      real*8 t(3,3,0:*),alpha,pi,angle
      character*(*) axis
c
      pi=4.0d+00*atan(1.0d+00)
      alpha=2*pi/naxis
c
      call rzero(t,3*3*naxis)
      if (axis.eq.' '.or.axis.eq.'z') then
         do 1 i=0,naxis-1
            angle=float(i)*alpha
            call rotz(t(1,1,i),angle)
    1    continue
      else if (axis.eq.'y') then
         do 2 i=0,naxis-1
            angle=float(i)*alpha
            call roty(t(1,1,i),angle)
    2    continue
      else if (axis.eq.'x') then
         do 3 i=0,naxis-1
            angle=float(i)*alpha
            call rotx(t(1,1,i),angle)
    3    continue
      else
         call lnkerr(' cnaxis: unrecognized axis')
      end if
c
c
      return
      end
