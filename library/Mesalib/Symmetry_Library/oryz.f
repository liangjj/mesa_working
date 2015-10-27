*deck @(#)oryz.f	5.1  11/6/94
      subroutine oryz(maxap3,a,b,natoms,atmchg,ixyz)
      implicit real*8(a-h,o-z)
c
c     the plane of a planar molecule is in the xy(ixyz=3) or
c     xz(ixyz=2) planes.  rotate it so that the molecular plane
c     is coincident with the yz cartesian plane.
c
      dimension a(maxap3,3), b(maxap3,3), atmchg(1), t(3,3), e(3)
      data one,two/1.0d+00,2.0d+00/
      save one,two
c
c     call rtrace(6horyz  ,1)
      numatm = natoms + 3
      halfpi = two * atan(one)
      goto(20,40,60), ixyz
c
 20   call orplan(maxap3,a,b,atmchg,numatm,e,t,1)
      return
c
 40   call rotate(maxap3,a,b,numatm,t,3,halfpi)
      call move(maxap3,b,a,numatm)
      call orplan(maxap3,a,b,atmchg,numatm,e,t,1)
      return
c
 60   call rotate(maxap3,a,b,numatm,t,2,halfpi)
      call move(maxap3,b,a,numatm)
      call orplan(maxap3,a,b,atmchg,numatm,e,t,1)
      return
      end
