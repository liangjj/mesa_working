*deck @(#)rcoord.f	5.1  11/6/94
      subroutine rcoord(file,natoms,ian,atchg,c)
c***begin prologue     rcoord
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           rwf, chk, coordinates,  atomic charges, atomic numbers
c***author             martin, richard (lanl)
c***source
c***purpose            retrieves the atomic numbers, charges, and coordinates
c                      from the internal file named 'file'.
c***description
c                      call rcoord(file,natoms,ian,atmchg,c)
c                        file     character string giving the file to search
c                                 (e.g., 'rwf').
c                        natoms   the number of atoms.
c                        ian      the atomic numbers.
c                        atmchg   the atomic charges.
c                        c        the atomic coordinates.
c
c***references
c
c***iosys i/o             unknown unit
c                      "atomic numbers"     integer      read
c                      "nuclear charges"     real         read
c                      coords        real         read
c
c***routines called    iosys(io)
c***end prologue       rcoord
      implicit integer(a-z)
      character*(*) file
      character filnam*8
      real*8 atchg(natoms),c(3*natoms)
      integer ian(natoms)
c
 
c
      filnam=file
      call iosys('read integer "atomic numbers" from '//filnam,
     $     -1,ian,0,' ')
      call iosys('read real "nuclear charges" from '//filnam,
     $     -1,atchg,0,' ')
      call iosys('read real coordinates from '//filnam,-1,c,0,' ')
c
c
      return
      end
