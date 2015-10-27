*deck @(#)tform.f	5.1  11/6/94
      subroutine tform(maxap3,t,a,b,n)
      implicit real*8(a-h,o-z)
c
c     t is the 3x3 transformation matrix which is used to transform
c     the n coordinates in a to those in b.
c
      dimension t(3,3), a(maxap3,3), b(maxap3,3)
c
c     call rtrace(6htform ,1)
      do 100 iat=1,n
         b(iat,1) = t(1,1)*a(iat,1) + t(1,2)*a(iat,2) + t(1,3)*a(iat,3)
         b(iat,2) = t(2,1)*a(iat,1) + t(2,2)*a(iat,2) + t(2,3)*a(iat,3)
         b(iat,3) = t(3,1)*a(iat,1) + t(3,2)*a(iat,2) + t(3,3)*a(iat,3)
 100  continue
      return
      end
