*deck @(#)wcoord.f	5.1  11/6/94
      subroutine wcoord(file,natoms,ian,atchg,c)
c***begin prologue     wcoord
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           rwf, chk, coordinates,  atomic charges, atomic numbers
c***author             martin, richard (lanl)
c***source
c***purpose            writes the atomic numbers, charges, and coordinates
c                      to the internal file named 'file'.
c***description
c                      call wcoord(file,natoms,ian,atmchg,c)
c                        file     character string giving the file to search
c                                 (e.g., 'rwf').
c                        natoms   the number of atoms.
c                        ian      the atomic numbers.
c                        atmchg   the atomic charges.
c                        c        the atomic coordinates.
c
c***references
c
c***iosys i/o                   unit unknown
c                      "atomic numbers"      integer     written
c                      "nuclear charges"      real        written
c                      coords         real        written
c***routines called    iosys(io)
c***end prologue       wcoord
      implicit integer(a-z)
      character*(*) file
      character filnam*8
      real*8 atchg(natoms),c(3*natoms)
      integer ian(natoms)
c
c
      filnam=file
      call iosys('write integer "atomic numbers" to '//filnam,
     $     natoms,ian,0,' ')
      call iosys('write real "nuclear charges" to '//filnam,
     $     natoms,atchg,0,' ')
      call iosys('write real coordinates to '//filnam,3*natoms,c,0,' ')
c
c
      return
      end
