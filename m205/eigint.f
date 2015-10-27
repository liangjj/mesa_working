*deck @(#)eigint.f	5.1  11/6/94
      subroutine eigint(eigvec,fre,natoms3,nvar,vname,b,evcint,
     $                  icfx)
c***begin prologue     eigint.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             page, michael(nrl)
c***source             @(#)eigint.f	5.1   11/6/94
c***purpose            transform cartesian vibrational eigenvectors
c                      to internals.
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       eigint.f
      implicit none
c     --- input variables -----
      integer nvar,natoms3,icfx
c     --- input arrays (unmodified) ---
      character*(*) vname(nvar)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 eigvec(natoms3,natoms3),evcint(nvar,nvar)
      real*8 fre(natoms3)
c     --- output variables ---
c     --- scratch arrays ---
      real*8 b(natoms3,nvar)
c     --- local variables ---
      integer inp,iout
      integer nrt,count,i,j,icol
      real*8 fac,zero,five
c
      parameter (zero=0.0d+00,five=5.0d+00)
      common /io/ inp,iout
c
 1000 format(5x,'cartesian normal modes ')
 1010 format(/
     $,' ** warning:  different $ of zero frequencies than expected'/)
 1020 format(5x,'internal coordinate normal modes ')
c
c     --- print the normal modes.
      write(iout,1000)
      call frqprt(eigvec,natoms3,natoms3,natoms3,natoms3,0,0,
     $            vname,' ',0,fre,.true.)
c
      if(icfx.eq.1) return
c
c     --- transform the eigenvectors to internal coordinates
c         first crunch the columns of the eigenvector matrix
      nrt=6
      count=0
      do 10 icol=1,natoms3
         if(abs(fre(icol)).gt.five) then
            count=count+1
            fre(count)=fre(icol)
            do 5 j=1,natoms3
               eigvec(j,count)=eigvec(j,icol)
    5       continue
         endif
   10 continue
      if(count.ne.natoms3-nrt) then
         write(iout,1010)
      endif
      call ebtc(evcint,b,eigvec,nvar,natoms3,nvar)
c
c     --- normalize the normal modes expressed in z-matrix coodinates
      do 200 j=1,nvar
         fac=zero
         do 100 i=1,nvar
            fac=fac+evcint(i,j)*evcint(i,j)
  100    continue
         do 150 i=1,nvar
            evcint(i,j)=evcint(i,j)/sqrt(fac)
  150    continue
  200 continue
c
c     --- print internal normal modes.
      write(iout,1020)
      call frqprt(evcint,nvar,nvar,nvar,nvar,1,0,vname,' ',
     $            0,fre,.true.)
c
c
      return
      end
