*deck @(#)basprt.f	1.1  11/20/92
      subroutine basprt(natoms,atsymb,atomz,basnam,typnam,nbtype,
     #                  ptprim,ptcont,nprim,ncont,ex,cf,numex,numcf,
     #                  start,iout,atomno,nctype,minmom,maxmom)
c
      implicit integer (a-z)
c
      character*(*) atsymb(0:*),basnam(natoms),typnam(nbtype)
      real*8 ex(numex),cf(numcf),atomz(natoms),pi32
      integer ptprim(natoms,nbtype),ptcont(natoms,nbtype)
      integer nprim(natoms,nbtype),ncont(natoms,nbtype)
      integer start(natoms,nbtype),atomno(natoms)
      integer nctype(nbtype),minmom(nbtype),maxmom(nbtype)
c
 1000 format(//,' basis set information.')
 1010 format(/' atom',i3,2x,a2/' charge:',f6.2/' basis:',a17)
 1020 format(' type:',a5,6x,'exponent',4x,'coefficient',
     $       i6,' primitives,',i5,' contracted functions.')
 1030 format(2x,a8,3(2x,f13.6))
 1050 format(' basis same as on atom',i4)
c
c     ----- start timing -----
c
c
c     ----- get pi**1.5 -----
c
      call iosys('read real pi from rwf',1,pi32,0,' ')
      pi32=pi32**1.5d+00
c
c
      write(iout,1000)
      do 100 atom=1,natoms
         write(iout,1010) atom,atsymb(atomno(atom)),atomz(atom),
     $                    basnam(atom)
         do 70 prvat=1,atom-1
            if(basnam(prvat).eq.basnam(atom)
     $         .and.atomno(atom).eq.atomno(prvat)) then
               write(iout,1050) prvat
               goto 100
            endif
   70    continue
         do 90 type=1,nbtype
            if(nprim(atom,type).gt.0) then
               write(iout,1020) typnam(type),nprim(atom,type),
     $                          ncont(atom,type)
               call baspr1(ex(ptprim(atom,type)),nprim(atom,type),
     $                     cf(ptcont(atom,type)),ncont(atom,type),
     $                     nctype(type),minmom(type),maxmom(type),
     $                     typnam,nbtype,iout,pi32)
            end if
   90    continue
  100 continue
c
c     ----- stop timing -----
c
c
c
      return
      end
