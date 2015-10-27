*deck @(#)normal.f	5.1  11/6/94
      subroutine normal(natoms,nbtype,ptprim,ptcont,nprim,
     $                  ncont,ex,numex,cf,numcf,iout,nctype,minmom,
     $                  maxmom)
c***begin prologue     normal
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)  .
c***keywords           one-electron, integrals, normalize
c***author             saxe, paul (lanl).
c***source
c***purpose            normalizes the contracted basis functions.
c***description
c                      call normal(natoms,nbtype,ptprim,ptcont,nprim,
c                                  ncont,ex,numex,cf,numcf,iout,nctype,
c                                  minmom,maxmom)
c***references
c***routines called    norm1(util)
c***end prologue       normal
      implicit integer (a-z)
c
      real*8 ex(numex),cf(numcf)
      integer ptprim(natoms,nbtype),ptcont(natoms,nbtype)
      integer nprim(natoms,nbtype),ncont(natoms,nbtype)
      integer nctype(nbtype),minmom(nbtype),maxmom(nbtype)
c
c     ----- start timing -----
c
c
c     ----- loop over atoms and types -----
c
      do 1000 atom=1,natoms
         do 900 type=1,nbtype
            if (nprim(atom,type).le.0) go to 900
c
c           ----- check if already normalized because of previous atom
c
c..bhl expanded if statement 8/8/89 llnl
            do 1 i=1,atom-1
               if (ptprim(i,type).eq.ptprim(atom,type) .and.
     $  nprim(i,type).ne.0) go to 1000
    1       continue
c
            call norm1(ex(ptprim(atom,type)),cf(ptcont(atom,type)),
     #                 nprim(atom,type),ncont(atom,type),nctype(type),
     #                 minmom(type),maxmom(type))
  900    continue
 1000 continue
c
c     ----- stop timing -----
c
c
c
      return
      end
