*deck @(#)zcentr.f	5.1  11/6/94
      subroutine zcentr(card,icur,outstr,ocur,iz,namcnt,nz)
c***begin prologue     zcentr.f
c***date written       850601  yymmdd
c***revision date      11/6/94
c***keywords           z-matrix
c***author             binkley, et.al., gaussian 82.
c***                   martin, richard (lanl)
c***source             @(#)zcentr.f	5.1   11/6/94
c***purpose            reads a center specification from a z-matrix card.
c***description
c     call zcentr(card,icur,outstr,ocur,iz,namcnt,nz)
c
c     module to read a center specification from a z-matrix card.
c     this specification may be either an integer (the sequential number
c     of a previous z-matrix card), or the name of a previously defined
c     center.
c
c        card     the z-matrix card being parsed.
c        icur     the parse begins at card(icur+1:)
c        outstr   an output string.
c        ocur     cursor for output string.
c        iz       the sequential number of the center being processed.
c        namcnt   a character array containing the names of
c                 the centers.
c        nz       the sequential number of the current z-matrix card.
c***references
c***routines called
c***end prologue       zcentr.f
      implicit none
c     --- input variables -----
      integer icur,iz,nz
c     --- input arrays (unmodified) ---
      character*(*) card,namcnt(*)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
      integer ocur
      character*(*) outstr
c     --- scratch arrays ---
c     --- local variables ---
      integer start,end
      integer lsubst,ctoi
      character found*16, ffnext*16, center*16
c
      found=ffnext(card,icur,start,end)
      center=card(start:end)
      if(found.eq.'string') then
         iz=lsubst(namcnt,center,nz-1)
      else if(found.eq.'integer') then
         iz=ctoi(center)
      else
         call lnkerr(' center specification must be an integer or a '
     $               //'string:'//center)
      endif
c
      outstr(ocur+1:)=namcnt(iz)(1:4)
      ocur=ocur+5
c
c
      return
      end
