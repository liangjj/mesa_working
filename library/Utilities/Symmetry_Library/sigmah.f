*deck @(#)sigmah.f	5.1  11/6/94
      subroutine sigmah(t,naxis,axis,told)
c
c     generate the transformation matrices for sigmah:
c     the default rotation axis is generally chosen as the z-axis,
c     in which case the sigmah is reflection through the xy-plane. if another
c     axis has been chosen, generate the appropriate cyclic permutation.
c
      implicit integer (a-z)
c
      real*8 t(3,3,naxis),told(3,3,naxis)
      character*(*) axis
c
      call vmove(t,told,3*3*naxis)
      if (axis.eq.' '.or.axis.eq.'z') then
c        sigmah=xy
         do 1 i=1,naxis
            t(3,3,i)=-t(3,3,i)
    1    continue
      else if (axis.eq.'y') then
c        sigmah=xz
         do 2 i=1,naxis
            t(2,2,i)=-t(2,2,i)
    2    continue
      else if (axis.eq.'x') then
c        sigmah=yz
         do 3 i=1,naxis
            t(1,1,i)=-t(1,1,i)
    3    continue
      else
         call lnkerr('sigmah: unrecognized principal axis')
      end if
c
c
      return
      end
