*deck @(#)vprod.f	5.1  11/6/94
      subroutine vprod(vp,x,y)
c***begin prologue     vprod
c***date written       850601  yymmdd
c***revision date      yymmdd  yymmdd
c***keywords           vector product, cross product
c***author             binkley,et.al., gauss82
c***source             @(#)vprod.f	5.1   11/6/94
c***purpose            computes the vector product: vp=x cross y .
c***description
c     call vprod(vp,x,y)
c       vp      the vector product(3).
c       x       the first vector(3).
c       y       the second vector(3).
c***references
c***routines called    (none)
c***end prologue       vprod
      implicit real*8(a-h,o-z)
c
      dimension vp(3),x(3),y(3)
c
      vp(1)=x(2)*y(3)-x(3)*y(2)
      vp(2)=x(3)*y(1)-x(1)*y(3)
      vp(3)=x(1)*y(2)-x(2)*y(1)
c
c
      return
      end
