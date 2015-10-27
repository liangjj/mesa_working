*deck @(#)symptr.f	3.1  11/17/92
      subroutine symptr(nprim,momatm,maxmom,natoms,nbtype,sympt,
     # mompt)
c
c***begin prologue     symptr
c***date written       880728   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           symmetry pointer
c***author             lengsfield, byron (llnl)
c***source             @(#)symptr.f	3.1   11/17/92
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       symptr
c
      implicit integer (a-z)
c
      integer nprim(natoms,nbtype)
      integer momatm(natoms)
      integer maxmom(nbtype)
      integer sympt(*),mompt(*)
c
      common/io/inp,iout
c
c     ----- associate a symmetry pointer with each basis function
c
      bf=0
      maxl=0
      do 20 atom=1,natoms
         write(iout,*)' atom ',atom
         do 10 type=1,nbtype
            npr=nprim(atom,type)
            if (npr.gt.0) then
             angmom=maxmom(type)
             write(iout,*)' type  bf angmom npr',type,bf,angmom,npr
             maxl=max(maxl,angmom)
              if    (angmom.eq.0) then
               do 11 i=1,npr
               bf=bf+1
               sympt(bf)=1
               mompt(bf)=0
  11           continue
              elseif(angmom.eq.1) then
               do 12 i=1,npr
               bf=bf+1
               sympt(bf)=2
               mompt(bf)=1
               bf=bf+1
               sympt(bf)=3
               mompt(bf)=1
               bf=bf+1
               sympt(bf)=1
               mompt(bf)=1
  12           continue
              elseif(angmom.eq.2) then
               do 13 i=1,npr
               bf=bf+1
               sympt(bf)=1
               mompt(bf)=2
               bf=bf+1
               sympt(bf)=1
               mompt(bf)=2
               bf=bf+1
               sympt(bf)=1
               mompt(bf)=2
               bf=bf+1
               sympt(bf)=4
               mompt(bf)=2
               bf=bf+1
               sympt(bf)=2
               mompt(bf)=2
               bf=bf+1
               sympt(bf)=3
               mompt(bf)=2
  13           continue
            elseif(angmom.eq.3) then
               do 14 i=1,npr
               bf=bf+1
               sympt(bf)=2
               mompt(bf)=3
               bf=bf+1
               sympt(bf)=3
               mompt(bf)=3
               bf=bf+1
               sympt(bf)=1
               mompt(bf)=3
               bf=bf+1
               sympt(bf)=2
               mompt(bf)=3
               bf=bf+1
               sympt(bf)=3
               mompt(bf)=3
               bf=bf+1
               sympt(bf)=1
               mompt(bf)=3
               bf=bf+1
               sympt(bf)=2
               mompt(bf)=3
               bf=bf+1
               sympt(bf)=3
               mompt(bf)=3
               bf=bf+1
               sympt(bf)=1
               mompt(bf)=3
               bf=bf+1
               sympt(bf)=4
               mompt(bf)=3
  14           continue
            else
               call lnkerr(' bug in symptr .. unknown maxmom ')
            end if
c
            end if
 10      continue
 20   continue
c
      write(iout,30) bf,maxl
  30  format(/,' m102: symptr has processed ',i6,' basis',
     # ' functions ',/,
     #         '       the maximum l-value was ',i6)
c
      call iosys('write integer sympt to rwf',bf,sympt,0,' ')
      call iosys('write integer mompt to rwf',bf,mompt,0,' ')
c
      return
      end
