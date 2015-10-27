*deck @(#)mulpop.f	5.1  11/6/94
      subroutine mulpop(nbf,natoms,s,d,totpop,atomp,atomc,orbpop,
     $                  ncont,nobf,nbtype,bfstrt,bfend,names,
     $                  ntodo,mintyp,nx,ny,nz,caltyp,ops,bflabl,refops,
     $                  nactiv)
c***begin prologue     mulpop.f
c***date written       860722  
c***revision date      11/6/94      
c
c***keywords           mulliken population, wave-function analysis
c***author             saxe, paul,    (lanl)
c***source             @(#)mulpop.f	5.1   11/6/94
c***purpose            the actual calculation and printing of hf
c                         mulliken populations.
c***description
c       mulpop forms and prints the mulliken populations for
c       closed- and open-shell rhf wavefunctions given the density matrices
c       and overlap matrix. for open-shell cases, the two density matrices
c       are the closed-shell and open-shell matrices. internally mulpop
c       calculates four copies of the population information for open-shell
c       cases: total, alpha, beta and difference populations, respectively.
c***references
c
c***routines called
c
c***end prologue       mulpop.f
      implicit none
c     --- input variables -----
      integer nbf,natoms,nbtype,ntodo,nactiv
      character*16 names(nactiv)
      character*(*) bflabl(*)
      character*(*) caltyp,ops
      character*(*) refops
c     --- input arrays (unmodified) ---
      integer ncont(natoms,nbtype),nobf(nbtype),bfstrt(nactiv)
      integer bfend(nactiv)
      integer mintyp(nbtype),nx(*),ny(*),nz(*)
c     --- input arrays (scratch) ---
      real*8 s(nbf,nbf),d(nbf,nbf,*),totpop(nbf,nbf,ntodo)
      real*8 atomc(nactiv,ntodo),orbpop(nbf,ntodo)
      real*8 atomp(nactiv,nactiv,ntodo)
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      character*16 collab(4)
      logical logkey
      integer inp,iout
      integer i,j,t,bf,atom,type,iatom,jatom
      integer ibf,jbf
      real*8 zero
c
      parameter (zero=0.0d+00)
c
      common /io/ inp,iout
c
    3 format(5x,'total mulliken population matrix:')
    6 format(5x,'gross orbital mulliken population:')
   26 format(5x,'spin-density matrix:')
   28 format(5x,'bond-order matrix:')
   31 format(5x,'total mulliken charges:',/,
     $        (10x,a8,2x,f11.5))
   32 format(5x,'total mulliken charges and spins:',/,
     $        (12x,a8,2f11.6))
   50 format(8x,a16,4(f9.6,3x))
c
c     --- form the total population matrix ---
      if (caltyp.eq.'rohf') then
         do 22 i=1,nbf
            do 21 j=1,nbf
               totpop(j,i,2)=s(j,i)*(d(j,i,1)+d(j,i,2))
               totpop(j,i,3)=s(j,i)*d(j,i,1)
               totpop(j,i,1)=totpop(j,i,2)+totpop(j,i,3)
               totpop(j,i,4)=totpop(j,i,2)-totpop(j,i,3)
   21       continue
   22    continue
      else if (caltyp.eq.'rhf') then
         do 2 i=1,nbf
            do 1 j=1,nbf
               totpop(j,i,1)=s(j,i)*d(j,i,1)*2.0d+00
    1       continue
    2    continue
      else if (caltyp.eq.'ci') then
         do 92 i=1,nbf
            do 91 j=1,nbf
               totpop(j,i,1)=s(j,i)*d(j,i,1)
   91       continue
   92    continue
      end if
      if(logkey(ops,'properties=populations=full',.false.,refops)) then
         write(iout,3)
         call wmat(totpop,nbf,nbf,bflabl,bflabl)
      end if
c
c     --- contract to orbitals ---
      do 23 t=1,ntodo
         do 5 i=1,nbf
            orbpop(i,t)=zero
            do 4 j=1,nbf
               orbpop(i,t)=orbpop(i,t)+totpop(j,i,t)
    4       continue
    5    continue
   23 continue
      if(logkey(ops,'properties=populations=orbital',.false.,refops))
     $     then
         write(iout,6)
         collab(1)='total'
         if(caltyp.eq.'rohf') then
            collab(2)='alpha'
            collab(3)='beta'
            collab(4)='spin'
         endif
         call wmat(orbpop,nbf,ntodo,bflabl,collab)
      endif
c
c     --- work out the starting and ending orbitals for each atom ---
      bf=0
      do 10 atom=1,nactiv
         bfstrt(atom)=bf+1
         do 9 type=1,nbtype
            bf=bf+ncont(atom,type)*nobf(type)
    9    continue
         bfend(atom)=bf
   10 continue
c
c     --- sum the total populations to atom centres ---
      do 25 t=1,ntodo
         do 20 iatom=1,nactiv
            do 19 jatom=1,nactiv
               atomp(iatom,jatom,t)=zero
               do 15 ibf=bfstrt(iatom),bfend(iatom)
                  do 14 jbf=bfstrt(jatom),bfend(jatom)
                     atomp(iatom,jatom,t)=atomp(iatom,jatom,t)+
     $                                    totpop(ibf,jbf,t)
   14             continue
   15          continue
   19       continue
   20    continue
   25 continue
      if(logkey(ops,'properties=populations=bond',.true.,refops)) then
         write(iout,28)
         call wmat(atomp,nactiv,nactiv,names,names)
         if(caltyp.eq.'rohf') then
            write(iout,26)
            call wmat(atomp(1,1,4),nactiv,nactiv,names,names)
         endif
      endif
c
c     --- sum the bond orders to get total charges ---
      do 27 t=1,ntodo
         do 30 iatom=1,nactiv
            atomc(iatom,t)=zero
            do 29 jatom=1,nactiv
               atomc(iatom,t)=atomc(iatom,t)+atomp(iatom,jatom,t)
   29       continue
   30    continue
   27 continue
      if(logkey(ops,'properties=populations=charge',.true.,refops))
     $     then
         if (caltyp.eq.'rohf') then
            write (iout,32) (names(i),atomc(i,1),
     $                      atomc(i,4),i=1,nactiv)
         else
            write (iout,31) (names(i),atomc(i,1),i=1,nactiv)
         end if
      end if
c
c
      return
      end
