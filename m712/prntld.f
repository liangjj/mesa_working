*deck @(#)prntld.f	5.1  11/6/94
      subroutine prntld(symcen,angmom,nprimi,nprimj,nprimk,npriml,
     $     nfi,nfj,nfk,nfl,istart,jstart,kstart,lstart,ndmat,
     $     dij,dik,dil,djk,djl,dkl)
c
c***begin prologue     prntld
c***date written       870703   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           local density matrices
c***author             saxe, paul (lanl)
c***source             @(#)prntld.f	5.1   11/6/94
c
c***purpose            printing local density matrices
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       prntld
c
      implicit integer (a-z)
c
      integer symcen(4)
      integer angmom(4)
      real*8 dij(nprimi,nprimj,nfi,nfj,ndmat)
      real*8 dik(nprimi,nprimk,nfi,nfk,ndmat)
      real*8 dil(nprimi,npriml,nfi,nfl,ndmat)
      real*8 djk(nprimj,nprimk,nfj,nfk,ndmat)
      real*8 djl(nprimj,npriml,nfj,nfl,ndmat)
      real*8 dkl(nprimk,npriml,nfk,nfl,ndmat)
c
      common /io/     inp,iout
c
c     ----- print out the density matrices, one-by-one -----
c
      write (iout,1000)
 1000 format (5x,'local density matrices:')
      write (iout,1010) symcen(1),angmom(1),symcen(2),
     $     angmom(2),nprimi,nprimj,nfi,nfj,istart,jstart
 1010 format (/,5x,'iatom imom jatom jmom nprimi nprimj ',
     $     'nfi nfj istart jstart',/,5x,i4,i5,i6,i5,i6,i7,
     $     i5,i4,i6,i7)
      do 1015 if=1,nfi
         do 1014 jf=1,nfj
            do 1013 dmat=1,ndmat
               write (iout,1011) if,jf,ndmat
 1011          format (5x,'if=',i2,'  jf=',i2,' shell=',i2)
               call matout(dij(1,1,if,jf,dmat),nprimi,nprimj,nprimi,
     $              nprimj,iout)
 1013       continue
 1014    continue
 1015 continue
c
      write (iout,1020) symcen(1),angmom(1),symcen(3),
     $     angmom(3),nprimi,nprimk,nfi,nfk,istart,kstart
 1020 format (/,5x,'iatom imom katom kmom nprimi nprimk ',
     $     'nfi nfk istart kstart',/,5x,i4,i5,i6,i5,i6,i7,
     $     i5,i4,i6,i7)
      do 1025 if=1,nfi
         do 1024 kf=1,nfk
            do 1023 dmat=1,ndmat
               write (iout,1021) if,kf,dmat
 1021          format (5x,'if=',i2,'  kf=',i2,' shell=',i2)
               call matout(dik(1,1,if,kf,dmat),nprimi,nprimk,nprimi,
     $              nprimk,iout)
 1023       continue
 1024    continue
 1025 continue
c
      write (iout,1030) symcen(1),angmom(1),symcen(4),
     $     angmom(4),nprimi,npriml,nfi,nfl,istart,lstart
 1030 format (/,5x,'iatom imom latom lmom nprimi npriml ',
     $     'nfi nfl istart lstart',/,5x,i4,i5,i6,i5,i6,i7,
     $     i5,i4,i6,i7)
      do 1035 if=1,nfi
         do 1034 lf=1,nfl
            do 1033 dmat=1,ndmat
               write (iout,1031) if,lf,dmat
 1031          format (5x,'if=',i2,'  lf=',i2,' shell=',i2)
               call matout(dil(1,1,if,lf,dmat),nprimi,npriml,nprimi,
     $              npriml,iout)
 1033       continue
 1034    continue
 1035 continue
c
      write (iout,1040) symcen(2),angmom(2),symcen(3),
     $     angmom(3),nprimj,nprimk,nfj,nfk,jstart,kstart
 1040 format (/,5x,'jatom jmom katom kmom nprimj nprimk ',
     $     'nfj nfk jstart kstart',/,5x,i4,i5,i6,i5,i6,i7,
     $     i5,i4,i6,i7)
      do 1045 jf=1,nfj
         do 1044 kf=1,nfk
            do 1043 dmat=1,ndmat
               write (iout,1041) jf,kf,dmat
 1041          format (5x,'jf=',i2,'  kf=',i2,' shell=',i2)
               call matout(djk(1,1,jf,kf,dmat),nprimj,nprimk,nprimj,
     $              nprimk,iout)
 1043       continue
 1044    continue
 1045 continue
c
      write (iout,1050) symcen(2),angmom(2),symcen(4),
     $     angmom(4),nprimj,npriml,nfj,nfl,jstart,lstart
 1050 format (/,5x,'jatom jmom latom lmom nprimj npriml ',
     $     'nfj nfl jstart lstart',/,5x,i4,i5,i6,i5,i6,i7,
     $     i5,i4,i6,i7)
      do 1055 jf=1,nfj
         do 1054 lf=1,nfl
            do 1053 dmat=1,ndmat
               write (iout,1051) jf,lf,dmat
 1051          format (5x,'jf=',i2,'  lf=',i2,' shell=',i2)
               call matout(djl(1,1,jf,lf,dmat),nprimj,npriml,nprimj,
     $              npriml,iout)
 1053       continue
 1054    continue
 1055 continue
c
      write (iout,1060) symcen(3),angmom(3),symcen(4),
     $     angmom(4),nprimk,npriml,nfk,nfl,kstart,lstart
 1060 format (/,5x,'katom kmom latom lmom nprimk npriml ',
     $     'nfk nfl kstart lstart',/,5x,i4,i5,i6,i5,i6,i7,
     $     i5,i4,i6,i7)
      do 1065 kf=1,nfk
         do 1064 lf=1,nfl
            do 1063 dmat=1,ndmat
               write (iout,1061) kf,lf,ndmat
 1061          format (5x,'kf=',i2,'  lf=',i2,' shell=',i2)
               call matout(dkl(1,1,kf,lf,dmat),nprimk,npriml,nprimk,
     $              npriml,iout)
 1063       continue
 1064    continue
 1065 continue
c
c
      return
      end
