*deck @(#)fcorpr.f	5.1  11/6/94
      subroutine fcorpr(iout,natoms,ian,f)
c***begin prologue     fcorpr.f
c***date written       850601  (yymmdd)
c***revision date      11/6/94
c***keywords           coordinates, forces, print
c***author             binkley, et al., gaussian82
c                      martin, richard (lanl)
c***source             @(#)fcorpr.f	5.1   11/6/94
c***purpose            prints cartesian coordinate forces.
c***description
c                      call fcorpr(iout,natoms,ian,f)
c                        iout     unit number of the output file.
c                        natoms   number of atoms.
c                        ian      integer atomic numbers,dimensioned (natoms).
c                        f        real array (3,natoms),
c                                 containing the (x,y,z) forces to print.
c
c***references
c***routines called    (none)
c***end prologue       fcorpr.f
      implicit none
c     --- input variables -----
      integer natoms,iout
c     --- input arrays (unmodified) ---
      integer ian(natoms)
      real*8 f(3,natoms)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer icent,iat
c
 1020 format(5x,'center     atomic           forces (hartrees/bohr)',
     $      /5x,'number     number         x           y           z')
 1030 format(5x,i3,7x,i4,4x,3f12.6)
 1040 format(5x,10x,  i4,4x,3f12.6)
c
c
      write(iout,1020)
      icent=0
      do 10 iat = 1, natoms
         if(ian(iat).ge.0) then
            icent=icent+1
            write(iout,1030) icent,ian(iat),f(1,iat),f(2,iat),f(3,iat)
         else
            write(iout,1040)       ian(iat),f(1,iat),f(2,iat),f(3,iat)
         endif
   10 continue
c
c
      return
      end
