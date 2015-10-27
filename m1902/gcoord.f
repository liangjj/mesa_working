*deck @(#)gcoord.f	5.1 11/6/94 
      subroutine gcoord(nat,c,inau,toang)
c***begin prologue     gcoord.f
c***date written       930222  
c***revision date      11/6/94      
c
c***keywords           cartesian, coordinates
c***author             martin, richard(lanl) 
c***source             @(#)gcoord.f	5.1   11/6/94
c***purpose            reads a list of cartesian coordinates
c                      from the input stream.
c***description
c    
c
c***references
c
c***routines called
c
c***end prologue       gcoord.f
      implicit none
c     --- input variables -----
      integer nat
      logical inau
      real*8 toang
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 c(*)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer iatom,cind,i,j,icur,bcur,ecur
      character*16 ffnext, found
      character*80 card
      real*8 ctofp
c
      common/io/inp,iout
c
 1000 format(a)
 1020 format(' looking for a cartesian coordinate and found: ',80a1)
c
c     --- module to read the geometry in cartesian coordinates.
c         loop over atoms until an end of section is found.
      iatom=0
      cind=0
   20 read(inp,1000,end=9000) card
         if(index(card,'$').ne.0) goto 40
         call locase(card,card)
         icur=0
         found=ffnext(card,icur,bcur,ecur)
c        --- skip blank cards.
         if(ecur.eq.0) goto 20
c
c        --- found something, now crack it, looking for three floating point 
c        entries.
         iatom=iatom+1
         icur=0
         do 30 i=1,3
            found=ffnext(card,icur,bcur,ecur)
            if(found.eq.'floating point') then
               cind=cind+1
               c(cind)=ctofp(card(bcur:ecur))
            else
               write(iout,1020) (card(j:j),j=bcur,ecur)
               call lnkerr(' ')
            endif
   30   continue
      goto 20
c
c     --- have found all coordinate input.
c         set nat, possibly convert to atomic units, and return.
c
   40 nat=iatom
      if(.not.inau) then
         do 50 i=1,cind
            c(i)=c(i)/toang
   50    continue
      endif
c
c
      return
 9000 call lnkerr(' no terminus found for the $grid section.')
      end
