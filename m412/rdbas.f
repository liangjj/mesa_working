*deck @(#)rdbas.f	5.1  11/6/94
      subroutine rdbas(filenm,ex,cont,c,ptprim,noprim,ptcont,nocont,
     $                 start,nctype,nocart,nobf,maxmom,minmom,mintyp,
     $                 nx,ny,nz,nat,nprim,ncont,ntypes,lenxyz)
c***begin prologue     rdbas
c***date written       850601  yymmdd
c***revision date      yymmdd  yymmdd
c***keywords           basis set, i/o, rwf, chk
c***author             martin, richard (lanl)
c***source             @(#)rdbas.f	5.1   11/6/94
c***purpose            retrieves basis set information from a specific file.
c***description
c     call rdbas(filenm,ex,cont,c,ptprim,noprim,ptcont,nocont,
c                start,nctype,nocart,nobf,maxmom,minmom,mintyp,
c                nxx,ny,nz,nat,nprim,ncont,ntypes,lenxyz)
c***references         (none)
c***routines called    (none)
c***end prologue       rdbas
      implicit integer(a-z)
      character*(*) filenm
      character file*16
      real*8 c(3,nat),ex(nprim),cont(ncont)
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(ntypes)
      integer nobf(ntypes),maxmom(ntypes),minmom(ntypes),mintyp(ntypes)
      integer nctype(ntypes),nx(lenxyz),ny(lenxyz),nz(lenxyz)
c
c
      file=filenm
      call iosys('read real exponents from '//file,-1,ex,0,' ')
      call iosys('read real "contraction coefficients" from '//
     $     file,-1,cont,0,' ')
      call iosys('read real coordinates from '//file,-1,c,0,' ')
      call iosys('read integer "pointer to primitives" from '//
     $     file,-1,ptprim,0,' ')
      call iosys('read integer "number of primitives" from '//
     $     file,-1,noprim,0,' ')
      call iosys('read integer "pointer to contraction coefficients"'//
     $     ' from '//file,-1,ptcont,0,' ')
      call iosys('read integer "number of contraction coefficients" '//
     $     'from '//file,-1,nocont,0,' ')
      call iosys('read integer "pointer to first function" from '//
     $     file,-1,start,0,' ')
      call iosys('read integer nctype from '//file,-1,nctype,0,' ')
      call iosys('read integer "number of cartesians" from '//
     $     file,-1,nocart,0,' ')
      call iosys('read integer "number of pure functions" from '//
     $     file,-1,nobf,0,' ')
      call iosys('read integer "minimum momentum" from '//
     $     file,-1,minmom,0,' ')
      call iosys('read integer "maximum momentum" from '//
     $     file,-1,maxmom,0,' ')
      call iosys('read integer "pointer to cartesians"from '//
     $     file,-1,mintyp,0,' ')
      call iosys('read integer "power of x" from '//file,-1,nx,0,' ')
      call iosys('read integer "power of y" from '//file,-1,ny,0,' ')
      call iosys('read integer "power of z" from '//file,-1,nz,0,' ')
c
c
      return
      end
