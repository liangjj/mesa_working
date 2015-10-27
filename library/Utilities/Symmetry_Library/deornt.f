*deck @(#)deornt.f	5.1  11/6/94
      subroutine deornt(rotmat,idir)
      implicit real*8(a-h,o-z)
c
c     given the rotation matrix rotmat, either rotate the symmetry
c     operations from the standard to the original orientation
c     (idir=0) or vice versa (idir=1).
c
      dimension rotmat(3,3), t(3,3)
c
      do 10 i = 1, 3
         do 10 j = 1, 3
            if(idir.eq.0) t(i,j) = rotmat(i,j)
            if(idir.ne.0) t(j,i) = rotmat(i,j)
 10      continue
c     call fixrep(t)
         return
         end
