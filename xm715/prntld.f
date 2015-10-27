*deck %W%  %G%
      subroutine prntld(symcen,angmom,nprimi,nprimj,nprimk,npriml,
     $                  nfi,nfj,nfk,nfl,istart,jstart,kstart,lstart,
     $                  ndmat,dij,dkl)
c***begin prologue     prntld.f
c***date written       870703   (yymmdd)
c***revision date      11/6/94
c
c   18 december, 1993  rlm at lanl
c      modifying to print only two matrices.
c***keywords           local density matrices
c***author             saxe, paul (lanl)
c***source             %W%   %G%
c
c***purpose            printing local density matrices
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       prntld.f
      implicit none
c     --- input variables -----
      integer nprimi,nprimj,nprimk,npriml
      integer nfi,nfj,nfk,nfl
      integer istart,jstart,kstart,lstart
      integer ndmat
c     --- input arrays (unmodified) ---
      integer symcen(4)
      integer angmom(4)
      real*8 dij(nprimi,nprimj,nfi,nfj,ndmat)
      real*8 dkl(nprimk,npriml,nfk,nfl,ndmat)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer if,jf,kf,lf,dmat
c
      common /io/ inp,iout
c
 1000 format (5x,'local density matrices:')
 1010 format (/,5x,'iatom imom jatom jmom nprimi nprimj ',
     $     'nfi nfj istart jstart',/,5x,i4,i5,i6,i5,i6,i7,
     $     i5,i4,i6,i7)
 1011 format (5x,'if=',i2,'  jf=',i2,' shell=',i2)
 1020 format (/,5x,'katom kmom latom lmom nprimk npriml ',
     $     'nfk nfl kstart lstart',/,5x,i4,i5,i6,i5,i6,i7,
     $     i5,i4,i6,i7)
 1021 format (5x,'kf=',i2,'  lf=',i2,' shell=',i2)
c
      write (iout,1000)
c
c     --- print out the density matrices, one-by-one ---
      write (iout,1010) symcen(1),angmom(1),symcen(2),
     $     angmom(2),nprimi,nprimj,nfi,nfj,istart,jstart
      do 50 if=1,nfi
         do 40 jf=1,nfj
            do 30 dmat=1,ndmat
               write (iout,1011) if,jf,ndmat
               call matout(dij(1,1,if,jf,dmat),nprimi,nprimj,nprimi,
     $              nprimj,iout)
   30       continue
   40    continue
   50 continue
c
c
      write (iout,1020) symcen(3),angmom(3),symcen(4),
     $     angmom(4),nprimk,npriml,nfk,nfl,kstart,lstart
      do 100 kf=1,nfk
         do 90 lf=1,nfl
            do 80 dmat=1,ndmat
               write (iout,1021) kf,lf,ndmat
               call matout(dkl(1,1,kf,lf,dmat),nprimk,npriml,nprimk,
     $              npriml,iout)
   80       continue
   90    continue
  100 continue
c
c
      return
      end
