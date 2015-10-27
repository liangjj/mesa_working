*deck @(#)corpr1.f	5.1  11/6/94
      subroutine corpr1(iout,natoms,ian,c,conver)
c***begin prologue     corpr1
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           coordinates, print
c***author             martin, richard (lanl)
c***source
c***purpose            coordinate printing routine.
c***description
c                      call corpr1(iout,natoms,ian,c,conver)
c                        iout     output unit number.
c                        natoms   number of atoms.
c                        ian      atomic numbers; integer(natoms).
c                        c        coordinates; real(3,natoms).
c                        conver   coordinates are converted to conver*c
c                                 before printing.
c
c***references
c***routines called    smul(math), fillel(util)
c***end prologue       corpr1
      implicit real*8(a-h,o-z)
      character el(106)*2
      real*8 c(3,natoms),cloc(3)
      integer ian(natoms)
c
      data maxel/104/, el(1)/'x'/
      save maxel,el
c
 1020 format(5x,'cd cent  el                   coordinates(angstroms)',
     $      /5x,'                            x           y           z'
     $      ,3f12.6)
 1030 format(5x,i2,2x,i2,3x,a2,9x,3f12.6)
 1040 format(5x,i2,4x,   3x,a2,9x,3f12.6)
c
c
c
      call fillel(0,maxel,el(2))
      write(iout,1020)
      icent=0
      do 10 iat=1,natoms
         call smul(cloc,c(1,iat),conver,3) 
         idx=ian(iat)+2
         if(ian(iat).gt.0) then
            icent=icent+1
            write(iout,1030) iat,icent,el(idx),
     $                       (cloc(i),i=1,3)
         else
            write(iout,1040) iat,      el(idx),
     $                       (cloc(i),i=1,3)
         endif
   10 continue
c
c
      return
      end
