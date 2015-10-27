*deck @(#)vec.f	5.1  11/6/94
      subroutine vec(small,ohoh,u,c,j,k)
c***begin prologue     vec
c***date written       850601  yymmdd
c***revision date      yymmdd  yymmdd
c***keywords           distance
c***author             gauss82
c***source             @(#)vec.f	5.1   11/6/94
c***purpose            computes length of vector connecting two centers.
c***description
c     call vec(small,ohoh,u,c,j,k)
c     module to compute the length of the vector connecting centers j,k.
c     if this is less than small, an error flag is set(ohoh).
c     if it is larger than small, the normalized vector is returned in u.
c***references
c***routines called    (none)
c***end prologue       vec
      implicit real*8(a-h,o-z)
      logical ohoh
      dimension c(1),r(3),u(3)
      data zero/0.0d+0/
      save zero
c
c
      r2=zero
      jtemp=(j-1)*3
      ktemp=(k-1)*3
      do 10 i=1,3
         r(i)=c(i+jtemp)-c(i+ktemp)
   10    r2=r2+r(i)*r(i)
      r2=sqrt(r2)
      ohoh = r2 .lt. small
      if (ohoh) return
      do 20 i=1,3
   20    u(i)=r(i)/r2
c
c
      return
      end
