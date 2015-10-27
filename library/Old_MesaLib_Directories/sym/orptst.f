*deck @(#)orptst.f	5.1  11/6/94
      subroutine orptst(maxap3,a,natoms,ixyz)
      implicit real*8(a-h,o-z)
c
c     if the molecule whose coordinates are in a is not planar or
c     is not contained in a cartesian plane, ixyz is set to zero.
c     otherwise it signifies the cartesian axis perpindicular to
c     the molecular plane.
c
      common/tol/ toler,tol2
      dimension a(maxap3,3)
c
c     call rtrace(6horptst,1)
      do 40 jxyz=1,3
         do 20 iat =1,natoms
            if (abs(a(iat,jxyz)) .gt. toler) goto 40
 20      continue
         ixyz = jxyz
         return
 40   continue
      ixyz = 0
      return
      end
