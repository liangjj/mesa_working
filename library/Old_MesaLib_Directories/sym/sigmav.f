*deck @(#)sigmav.f	5.1  11/6/94
      subroutine sigmav(t,naxis,axis,told)
c
c     generate the sigmav transformation matrices:
c     the default principal axis is the z-axis, which implies reflection
c     through the x-z plane. if another principal axis was chosen, generate
c     sigmav from cyclic permutations.
c
      implicit integer (a-z)
c
      real*8 t(3,3,naxis),told(3,3,naxis)
      character*(*) axis
c
      call vmove(t,told,3*3*naxis)
c
      if (axis.eq.' '.or.axis.eq.'z') then
c        sigmav=xz
         do 2 i=1,naxis
            do 1 j=1,3
               t(2,j,i)=-t(2,j,i)
    1       continue
    2    continue
      else if (axis.eq.'y') then
c        sigmav=yz
         do 4 i=1,naxis
            do 3 j=1,3
               t(1,j,i)=-t(1,j,i)
    3       continue
    4    continue
      else if (axis.eq.'x') then
c        sigmav=xy
         do 6 i=1,naxis
            do 5 j=1,3
               t(3,j,i)=-t(3,j,i)
    5       continue
    6    continue
      else
         call lnkerr('sigmav: unrecognized axis')
      end if

c
c
      return
      end
