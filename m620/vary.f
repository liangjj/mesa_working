*deck @(#)vary.f	5.1 11/6/94
      subroutine vary(ops,nbond,nat,ndip,nsave,nsite,nq,nbegin,
     $                jz,site,zl,tsite,ident,label,isel)
c***begin prologue     vary.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             williams, d.e.
c***source             @(#)vary.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       vary.f
      implicit none
c     --- input variables -----
      character*(*) ops
      character*(*) label(*)
      integer nbond,nat,nsave
      integer jz(nat),ident(2,nat)
c     --- input arrays (unmodified) ---
      real*8 site(3,nat)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      integer isel(nat,10)
      integer ndip(5,nsave)
      real*8 zl(3,nat)
c     --- output variables ---
      integer nsite,nq,nbegin
c     --- scratch arrays ---
      real*8 tsite(3,nat)
c     --- local variables ---
      logical logkey,monopole,dipoles,quads
      integer i,j,k,l
      integer inp,iout
c
      common/io/inp,iout
c
 1000 format(1x,5i5)
c
      nsite=nat
      if (nbond.eq.0) then
c        --- set up a mask which says whether to use charges,dipoles,etc.
c            this will not be as flexible as the original code, but can
c            be easily changed if need be.
         monopole=logkey(ops,'cox-williams=charges',.true.,' ')
         dipoles=logkey(ops,'cox-williams=dipoles',.true.,' ')
         quads=logkey(ops,'cox-williams=quadrupoles',.true.,' ')
         do 30 i=1,nat
            if(monopole) then
               isel(i,1)=1
            endif
            if(dipoles) then
               isel(i,2)=1
               isel(i,3)=1
               isel(i,4)=1
            endif
            if(quads) then
               isel(i,5)=0
               isel(i,6)=1
               isel(i,7)=1
               isel(i,8)=1
               isel(i,9)=1
               isel(i,10)=1
            endif
   30    continue
      else if(nbond.eq.1) then
c        --- set up bond dipole calculation; note that nsite is reset in
c            subroutine bond
         nbegin=1
         call bond(jz,site,zl,nsite,ident,nbegin,tsite)
c        --- set isel for the three bond dipole components
         do 60 j=1,nsite
            isel(j,1)=0
            do 40 i=2,4
               isel(j,i)=1
   40       continue
            do 50 i=5,10
               isel(j,i)=0
   50       continue
   60    continue
      else if (nbond.eq.2) then
c        --- set up restricted bond dipole calculation
         nbegin=1
         if (nsave.ne.0) then
            call intarr(ops,'cox-williams=ndip',ndip,5*nsave,' ')
            write(iout,1000) ((ndip(i,j),i=1,5),j=1,nsave)
c           --- insert atomic dipoles at beginning of list 
c               and find directions
            do 100 j=1,nsave
c              --- move everything up one place
               do 80 k=nsite,1,-1
                  jz(nbegin+k)=jz(nbegin+k-1)
                  label(nbegin+k)=label(nbegin+k-1)
                  do 70 i=1,3
                     site(i,nbegin+k)=site(i,nbegin+k-1)
   70             continue
   80          continue
               nbegin=nbegin+1
c              --- move the selected atom l into position j
               l=ndip(2,j)+j
               jz(j)=jz(l)
               label(j)=label(l)
               ident(1,j)=j
               ident(2,j)=j
               do 90 i=1,3
                  site(i,j)=site(i,l)
   90          continue
c                 now find the dipole direction zl of atom j
c                 note that in site(i,l) l has been increased by j
               call lonep(ndip,nsite,site,zl,j)
  100       continue
         endif
         call bond(jz,site,zl,nsite,ident,nbegin,tsite)
c        --- set isel for one dipole magnitude
c**** this comment doesn't sem to go with the code.
         do 120 j=1,nsite
            isel(j,1)=1
            do 110 i=2,10
               isel(j,i)=0
  110       continue
  120    continue
      endif
c
c     --- count variables
      nq=0
      do 300 l=1,10
         do 290 i=1,nsite
            if (isel(i,l).eq.1) nq=nq+1
c           --- test for bad input 
            if (isel(i,5).ne.0) then
               call lnkerr('pm620: do not vary qxx')
            endif
  290    continue
  300 continue
c
c
      return
      end
