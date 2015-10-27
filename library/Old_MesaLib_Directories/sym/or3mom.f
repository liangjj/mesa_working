*deck @(#)or3mom.f	5.1  11/6/94
      function or3mom(maxap3,a,atmchg,natoms,ixyz)
      implicit real*8(a-h,o-z)
c
c     this function returns the value of the third moment of charge
c     along the ixyz cartesian axis.  note that the distance used in
c     computing the moment is not as usually defined.  rather than
c     being the perpindicular distance from the axis to the point it
c     is the projection of the point onto the reference axis.
c
      dimension a(maxap3,3), atmchg(1)
      data zero/0.0d+00/
      save zero
c
c     call rtrace(6hor3mom,1)
      or3mom = zero
      do 100 iat=1,natoms
 100     or3mom = or3mom + atmchg(iat)*a(iat,ixyz)**3
         return
         end
