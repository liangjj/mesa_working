*deck ssssig
      subroutine ssssig(maxap3,natoms,iatflg,a,ixyz,itst)
      implicit real*8(a-h,o-z)
c
c     flag those atoms which are in the plane of symmetry perpindicular
c     to cartesian axes ixyz at ixyz=0.
c
      common/tol/ toler,tol2
      dimension iatflg(1), a(maxap3,3)
c
c     call rtrace(6hssssig,1)
      itst = 0
      do 20 iat=1,natoms
         if (abs(a(iat,ixyz)) .gt. toler  .or.
     $        iatflg(iat) .ne. 0) goto 20
         iatflg(iat) = 2
         itst = 1
 20   continue
      return
      end
