*deck @(#)ecpprt.f	1.1  11/20/92
      subroutine ecpprt(natoms,atsymb,atomz,basnam,typnam,ntypes,nbtype,
     #                  ptprim,ptcont,nprim,ncont,ex,cf,numex,numcf,
     #                  start,iout,atomno,nctype,minmom,maxmom)
c
      implicit integer (a-z)
c
      character*(*) atsymb(0:*),basnam(natoms),typnam(ntypes)
      real*8 ex(numex),cf(numcf),atomz(natoms)
      integer ptprim(natoms,ntypes),ptcont(natoms,ntypes)
      integer nprim(natoms,ntypes),ncont(natoms,ntypes)
      integer start(natoms,ntypes),atomno(natoms)
      integer nctype(ntypes),minmom(ntypes),maxmom(ntypes)
c
 1000 format(//,' effective core potential information.')
 1010 format(/' atom',i3,2x,a2/' charge:',f6.2/' basis:',a17)
 1020 format(' projector       exponent    coefficient           r**n')
 1030 format(2x,a8,3(2x,f16.6))
 1040 format(' no core electrons replaced.')
 1050 format(' ecp same as on atom',i4)
c
c     start timing.
c
      write(iout,1000)
      do 100 atom=1,natoms
         write(iout,1010) atom,atsymb(atomno(atom)),atomz(atom),
     $                    basnam(atom)
         if(atomno(atom).eq.int(atomz(atom))) then
            write(iout,1040)
            goto 100
         endif
         do 70 prvat=1,atom-1
            if(atomno(prvat).eq.atomno(atom)) then
               if(basnam(prvat).eq.basnam(atom)) then
                  write(iout,1050) prvat
                  goto 100
               endif
            endif
   70    continue
         write(iout,1020)
         do 90 type=nbtype+1,ntypes
            if(nprim(atom,type).gt.0) then
               do 80 j=0,nprim(atom,type)-1
                  write(iout,1030) typnam(type),ex(ptprim(atom,type)+j),
     $                         cf(ptcont(atom,type)+j),
     $                         cf(ptcont(atom,type)+nprim(atom,type)+j)
   80          continue
            end if
   90    continue
  100 continue
c
c     stop timing.
c
c
      return
      end
