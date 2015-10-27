*deck @(#)cprint.f	5.1  11/6/94
      subroutine cprint(natoms,ian,zan,c,conver,bastyp)
c***begin prologue     cprint.f
c***date written       850601  (yymmdd)
c***revision date      11/6/94
c***keywords           coordinates, print
c***author             martin, richard (lanl)
c***source             @(#)cprint.f	5.1   11/6/94
c***purpose            coordinate printing routine.
c***description
c                      call cprint(natoms,ian,z,c,conver,bastyp)
c                        natoms   number of atoms.
c                        ian      atomic numbers; integer(natoms).
c                        zan      atomic charges; real(natoms).
c                        c        coordinates; real(3,natoms).
c                        conver   coordinates are converted to conver*c
c                                 before printing.
c                        bastyp   basis set array; character*16 (natoms)
c
c***references
c***routines called    smul(math), fillel(util)
c***end prologue       cprint.f
      implicit none
c     --- input variables -----
      integer natoms
      real*8 conver
c     --- input arrays (unmodified) ---
      integer ian(natoms)
      character*(*) bastyp(natoms)
      real*8 c(3,natoms),zan(natoms)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer maxel
      integer icent,iat,i,idx
      character el(106)*2
      real*8 cloc(3)
c
      data maxel/104/, el(1)/'x'/
      save maxel,el
c
      common/io/inp,iout
c
 1010 format(1x,'cartesian coordinates(angstroms):')
 1020 format(5x,' cd  cent  el basis          z',
     $       t44,'coordinates(angstroms)',
     $      /13x,'                            x           y           z'
     $      )
 1030 format(5x,i3,2x,i3,3x,a2,1x,a8,f8.3,3f12.6)
 1040 format(5x,i3,5x,   3x,a2,1x,a8,f8.3,3f12.6)
c
c     --- fill the array of element names.
      call fillel(0,maxel,el(2))
      write(iout,1010)
      write(iout,1020)
      icent=0
      do 10 iat=1,natoms
         call smul(cloc,c(1,iat),conver,3)
         idx=ian(iat)+2
         if(ian(iat).gt.0) then
            icent=icent+1
            write(iout,1030) iat,icent,el(idx),bastyp(iat),zan(iat),
     $                       (cloc(i),i=1,3)
         else
            write(iout,1040) iat,      el(idx),bastyp(iat),zan(iat),
     $                       (cloc(i),i=1,3)
         endif
   10 continue
c
c
      return
      end
